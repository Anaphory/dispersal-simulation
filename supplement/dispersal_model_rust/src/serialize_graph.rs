use model::ecology;
use model::movementgraph::MovementGraph;
use std::collections::HashMap;
use std::fs::File;
use std::io::prelude::*;

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

    let mut graph: MovementGraph = MovementGraph::default();

    let mut nodes_stmt;
    match conn.prepare(concat!(
        "SELECT hexbin, longitude, latitude, popdensity FROM hex ",
        "WHERE ? < longitude AND longitude < ? AND ? < latitude AND latitude < ? ",
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
    match conn.prepare("SELECT ecoregion, area_m2 FROM eco WHERE hexbin = ? AND ecoregion != 999") {
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
        "WHERE ? < longitude AND longitude < ? AND ? < latitude AND latitude < ? ",
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

    /*
    let data = match ecology::load_density_tif() {
        Some(d) => d,
        None => return Err("Population density raster not found".to_string())
    };
    */
    let mut h3_to_graph = HashMap::new();
    let mut expand_attested = HashMap::new();
    println!("Adding nodes…");
    for (i, (hexbin, longitude, latitude, popdensity)) in nodes_stmt
        .query_map(
            rusqlite::params![boundary_west, boundary_east, boundary_south, boundary_north],
            |node| {
                Ok((
                    node.get::<_, i64>(0)?,
                    node.get(1)?,
                    node.get(2)?,
                    node.get::<_, f64>(3)?,
                ))
            },
        )
        .unwrap()
        .flatten()
        .enumerate()
    {
        println!("{}: {}", i, hexbin);
        let ecos: HashMap<_, _> = eco_stmt
            .query_map(rusqlite::params![hexbin], |eco| {
                // The areas are already scaled by the cosine of latitude,
                // so this converts them into km².
                let len = expand_attested.len();
                let ecoregion = eco.get::<_, i64>(0)?;
                Ok((
                    *(expand_attested.entry(ecoregion).or_insert(len)),
                    // Population density is in individuals per 100km², area is in m²
                    popdensity * eco.get::<_, f64>(1)? / 100_000_000.,
                ))
            })
            .unwrap()
            .flatten()
            .collect();
        // assert!(!ecos.is_empty());
        let h = hexbin;
        h3_to_graph.insert(
            h,
            graph.add_node((h as u64, longitude, latitude, ecos)),
        );
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
        .filter_map(|(i, j, d)| Some((*h3_to_graph.get(&i)?, *h3_to_graph.get(&j)?, d)))
    {
        graph.add_edge(a1, a2, d);
    }
    println!("{:} edges added.", graph.edge_count());
    Ok(graph)
}

fn parse_args(west: &mut f64, east: &mut f64, south: &mut f64, north: &mut f64) {
    let mut parser = argparse::ArgumentParser::new();
    parser.set_description("Run a dispersal simulation");
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

    parser.parse_args_or_exit();
}

fn main() -> Result<(), String> {
    let mut boundary_west = -168.571_541;
    let mut boundary_east = -34.535_395;
    let mut boundary_south = -56.028_198;
    let mut boundary_north = 74.52671;

    parse_args(
        &mut boundary_west,
        &mut boundary_east,
        &mut boundary_south,
        &mut boundary_north,
    );

    let dispersal_graph =
        read_graph_from_db(boundary_west, boundary_east, boundary_south, boundary_north)?;

    let mut file = match File::create("graph.bincode") {
        Ok(f) => f,
        Err(e) => return Err(e.to_string()),
    };
    match file.write(&bincode::serialize(&dispersal_graph).unwrap()) {
        Ok(_) => Ok(()),
        Err(e) => Err(e.to_string()),
    }
}
