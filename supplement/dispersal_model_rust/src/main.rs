use model::hexgrid;
use model::movementgraph::MovementGraph;
use model::submodels::parameters::Parameters;
use model::*;
use std::collections::HashMap;

fn read_graph_from_db(
    boundary_west: f64,
    boundary_east: f64,
    boundary_south: f64,
    boundary_north: f64,
) -> Result<MovementGraph, String> {
    let try_conn = rusqlite::Connection::open_with_flags(
        "/home/gereon/Public/settlement-of-americas/supplement/distances/plot.sqlite",
        rusqlite::OpenFlags::SQLITE_OPEN_READ_ONLY,
    );
    let conn: rusqlite::Connection;
    match try_conn {
        Err(e) => {
            return Err(e.to_string());
        }
        Ok(k) => {
            conn = k;
        }
    }

    let mut graph: MovementGraph = petgraph::csr::Csr::new();

    let mut nodes_stmt;
    match conn.prepare(concat!(
        "SELECT hexbin, vlongitude, vlatitude FROM hex ",
        "WHERE ? < vlongitude AND vlongitude < ? AND ? < vlatitude AND vlatitude < ? ",
        // "NATURAL INNER JOIN (SELECT hexbin FROM eco WHERE ecoregion != 999 GROUP BY hexbin) ",
    )) {
        Err(e) => {
            return Err(e.to_string());
        }
        Ok(k) => {
            nodes_stmt = k;
        }
    }

    let mut eco_stmt;
    match conn.prepare("SELECT ecoregion, frequency FROM eco WHERE hexbin = ? AND ecoregion != 999")
    {
        Err(e) => {
            return Err(e.to_string());
        }
        Ok(k) => {
            eco_stmt = k;
        }
    }

    let mut dist_stmt;
    match conn.prepare(concat!(
        "SELECT hexbin1, hexbin2, min(distance) ",
        "FROM dist ",
        "JOIN hex ON hexbin = hexbin1 ",
        "WHERE ? < vlongitude AND vlongitude < ? AND ? < vlatitude AND vlatitude < ? ",
        // "INNER JOIN (SELECT hexbin FROM eco WHERE ecoregion != 999 GROUP BY hexbin) ON hexbin = hexbin1",
        // "INNER JOIN (SELECT hexbin FROM eco WHERE ecoregion != 999 GROUP BY hexbin) ON hexbin = hexbin2",
        "GROUP BY hexbin1, hexbin2 ORDER BY hexbin1, hexbin2 "
    )) {
        Err(_) => {
            return Err("Could not prepare dist statement".to_string());
        }
        Ok(k) => {
            dist_stmt = k;
        }
    }

    let mut h3_to_graph = HashMap::new();
    let mut expand_attested = HashMap::new();
    println!("Adding nodes…");
    for (i, (hexbin, longitude, latitude)) in nodes_stmt
        .query_map(
            rusqlite::params![boundary_west, boundary_east, boundary_south, boundary_north],
            |node| Ok((node.get::<_, i64>(0)?, node.get(1)?, node.get(2)?)),
        )
        .unwrap()
        .flatten()
        .enumerate()
    {
        let ecos: HashMap<_, _> = eco_stmt
            .query_map(rusqlite::params![hexbin], |eco| {
                // The areas are already scaled by the cosine of latitude,
                // so this converts them into km².
                let len = expand_attested.len();
                Ok((
                    *(expand_attested.entry(eco.get::<_, i64>(0)?).or_insert(len)),
                    eco.get::<_, f64>(1)? as f32,
                ))
            })
            .unwrap()
            .flatten()
            .collect();
        // assert!(!ecos.is_empty());
        graph.add_node((hexbin as hexgrid::Index, longitude, latitude, ecos));
        h3_to_graph.insert(hexbin, i);
    }
    println!("{:} nodes added.", graph.node_count());
    assert!(expand_attested.len() <= ecology::ATTESTED_ECOREGIONS);
    println!("{:} ecoregions found.", expand_attested.len());
    println!("{:?}", expand_attested);

    println!("Adding edges…");
    for (a1, a2, d) in dist_stmt
        .query_map(
            rusqlite::params![boundary_west, boundary_east, boundary_south, boundary_north],
            |dist| Ok((dist.get(0)?, dist.get(1)?, dist.get(2)?)),
        )
        .unwrap()
        .flatten()
        .filter_map(|(i, j, d)| {
            print!("{:}", i);
            Some((*h3_to_graph.get(&i)?, *h3_to_graph.get(&j)?, d))
        })
    {
        graph.add_edge(a1, a2, d);
    }
    println!("{:} edges added.", graph.edge_count());
    Ok(graph)
}

fn main() -> Result<(), String> {
    let mut boundary_west = -168.571_541;
    let mut boundary_east = -34.535_395;
    let mut boundary_south = -56.028_198;
    let mut boundary_north = 74.52671;

    let mut p = Parameters {
        attention_probability: 0.1,
        time_step_energy_use: 1., // 2263. as KCal * 365.242_2 / 2.,
        storage_loss: 0.33,
        resource_recovery: 0.20,
        culture_mutation_rate: 6e-3,
        culture_dimensionality: 20,
        cooperation_threshold: 6,
        cooperation_gain: 0.5,
        accessible_resources: 0.2,
        evidence_needed: 0.3,
        payoff_std: 0.1,

        dispersal_graph: MovementGraph::new(),
    };
    let mut max_t: HalfYears = 20000;
    parse_args(
        &mut boundary_west,
        &mut boundary_east,
        &mut boundary_south,
        &mut boundary_north,
        &mut p,
        &mut max_t,
    );
    println!("# Initialization ...");
    // FIXME: There is something messed up in lat/long -> image pixels. It works
    // currently, but the names point out that I misunderstood something, eg.
    // with the layout of the tiff in memory or similar.

    p.dispersal_graph =
        read_graph_from_db(boundary_west, boundary_east, boundary_south, boundary_north)?;
    let s: State = initialization(&p).unwrap();
    println!("Initialized");

    run(s, p, max_t);
    Ok(())
}

fn parse_args(
    west: &mut f64,
    east: &mut f64,
    south: &mut f64,
    north: &mut f64,
    p: &mut Parameters,
    max_t: &mut HalfYears,
) {
    let mut parser = argparse::ArgumentParser::new();
    parser.set_description("Run a dispersal simulation");
    // TODO: Attach an ArgumentParser to the parameters object, so parameters can be set from the CLI.
    parser.refer(&mut p.attention_probability).add_option(
        &["--attention-probability"],
        argparse::Store,
        "attention probability",
    );
    parser.refer(&mut p.culture_mutation_rate).add_option(
        &["--culture-mutation-rate"],
        argparse::Store,
        "culture mutation rate",
    );
    parser.refer(&mut p.culture_dimensionality).add_option(
        &["--culture-dimensionality"],
        argparse::Store,
        "culture dimensionality",
    );
    parser.refer(&mut p.cooperation_threshold).add_option(
        &["--cooperation-threshold"],
        argparse::Store,
        "threshold under which cooperation happens",
    );
    parser.refer(&mut p.cooperation_gain).add_option(
        &["--cooperation-gain"],
        argparse::Store,
        "exponent in benefits from cooperation",
    );
    parser.refer(max_t).add_option(
        &["--steps"],
        argparse::Store,
        "number of half years to simulate",
    );
    parser.refer(west).add_option(
        &["--west"],
        argparse::Store,
        "western bounding box boundary",
    );
    parser.refer(east).add_option(
        &["--east"],
        argparse::Store,
        "eastern bounding box boundary",
    );
    parser.refer(south).add_option(
        &["--south"],
        argparse::Store,
        "southern bounding box boundary",
    );
    parser.refer(north).add_option(
        &["--north"],
        argparse::Store,
        "northern bounding box boundary",
    );
    println!("Starting…");
    parser.parse_args_or_exit();
}
