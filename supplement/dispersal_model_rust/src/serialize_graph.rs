use itertools::Itertools;
use model::ecology;
use model::movementgraph::MovementGraph;
use petgraph::visit::EdgeRef;
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
        "/home/gereon/Public/settlement-of-americas/supplement/compute_distances/network.sqlite",
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
        "SELECT node_id, longitude, latitude FROM nodes ",
        "WHERE ? < longitude AND longitude < ? AND ? < latitude AND latitude < ? ",
    )) {
        Err(e) => {
            return Err(e.to_string());
        }
        Ok(k) => {
            nodes_stmt = k;
        }
    }

    let mut eco_stmt;
    match conn.prepare("SELECT ecoregion, area, population_capacity FROM ecology WHERE node = ? AND ecoregion != 999") {
        Err(e) => {
            return Err(e.to_string());
        }
        Ok(k) => {
            eco_stmt = k;
        }
    }

    let mut dist_stmt;
    match conn.prepare(concat!(
        "SELECT node1, node2, min(travel_time) ",
        "FROM edges ",
        "JOIN nodes ON node1 = node_id ",
        "WHERE ? < longitude AND longitude < ? AND ? < latitude AND latitude < ? ",
        "GROUP BY node1, node2 ORDER BY node1, node2 "
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
    for (i, (node_id, longitude, latitude)) in nodes_stmt
        .query_map(
            rusqlite::params![boundary_west, boundary_east, boundary_south, boundary_north],
            |node| Ok((node.get::<_, i64>(0)?, node.get(1)?, node.get(2)?)),
        )
        .unwrap()
        .flatten()
        .enumerate()
    {
        println!("{}: {}", i, node_id);
        let ecos: HashMap<_, _> = eco_stmt
            .query_map(rusqlite::params![node_id], |eco| {
                let ecoregion = eco.get::<_, i64>(0)?;
                let index_if_missing = expand_attested.len();
                Ok((
                    *(expand_attested.entry(ecoregion).or_insert(index_if_missing)),
                    eco.get::<_, f64>(1)?,
                ))
            })
            .unwrap()
            .flatten()
            .collect();
        // assert!(!ecos.is_empty());
        h3_to_graph.insert(
            node_id,
            graph.add_node((node_id as u64, longitude, latitude, ecos)),
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

    let mut dispersal_graph =
        read_graph_from_db(boundary_west, boundary_east, boundary_south, boundary_north)?;

    for node in dispersal_graph.node_indices() {
        let node_attribs = dispersal_graph.node_weight(node);
        let population = match node_attribs {
            None => continue,
            Some(n) => (n.3.values().sum::<f64>()),
        };
        if population > 0.01 {
            println!("Population of node {:?}: {}", node, population);
            continue;
        }
        let new_edges: Vec<_> = dispersal_graph
            .edges_directed(node, petgraph::Direction::Incoming)
            .cartesian_product(dispersal_graph.edges_directed(node, petgraph::Direction::Outgoing))
            .filter_map(|(in_edge, out_edge)| {
                let weight = in_edge.weight() + out_edge.weight();
                if (in_edge.source() != out_edge.target()) & (weight < 3. * 8. * 3600.) {
                    Some((in_edge.source(), out_edge.target(), weight))
                } else {
                    None
                }
            })
            .collect();
        for (source, target, weight) in new_edges {
            match dispersal_graph.find_edge(source, target) {
                None => {
                    dispersal_graph.add_edge(source, target, weight);
                }
                Some(e) => {
                    if dispersal_graph.edge_weight(e).unwrap() > &weight {
                        dispersal_graph.update_edge(source, target, weight);
                    }
                }
            }
        }
        dispersal_graph.remove_node(node);
        println!("Reduced node {:?}", node);
    }

    let mut file = match File::create("graph.bincode") {
        Ok(f) => f,
        Err(e) => return Err(e.to_string()),
    };
    match file.write(&bincode::serialize(&dispersal_graph).unwrap()) {
        Ok(_) => Ok(()),
        Err(e) => Err(e.to_string()),
    }
}
