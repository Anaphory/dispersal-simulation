use model::argparse::parse_args;
use model::hexgrid;
use model::movementgraph::MovementGraph;
use model::submodels::parameters::Parameters;
use model::*;
use std::collections::HashMap;
use std::fs;

fn main() -> Result<(), String> {
    let contents = match fs::read("graph.bincode") {
        Ok(c) => c,
        Err(e) => return Err(e.to_string()),
    };

    let dispersal_graph: MovementGraph = match bincode::deserialize(&contents) {
        Ok(c) => c,
        Err(e) => return Err(e.to_string()),
    };

    let mut p = Parameters::default();
    p.dispersal_graph = dispersal_graph;
    let mut max_t: HalfYears = 20000;
    let mut scale: KCal = 1.0;
    parse_args(&mut p, &mut max_t, &mut scale);
    println!("# Initialization ...");
    // FIXME: There is something messed up in lat/long -> image pixels. It works
    // currently, but the names point out that I misunderstood something, eg.
    // with the layout of the tiff in memory or similar.

    let s: State = initialization(&p, &scale).unwrap();
    println!("Initialized");

    run(s, p, max_t);
    Ok(())
}
