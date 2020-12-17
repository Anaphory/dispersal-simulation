use model::argparse::parse_args;
use model::movementgraph::MovementGraph;
use model::submodels::parameters::Parameters;
use model::*;
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
    let mut max_t: Seasons = 20000;
    let mut scale = 1.0;
    let mut recovery = p.resource_recovery_per_season / p.season_length_in_years;

    let mut o = observation::ObservationSettings {
        log_every: 1,
        log_gdcd: 0,
        log_patch_resources: 1,
        store_every: 100,
        statefile: "state.json".to_string(),
    };

    {
        let parser = parse_args(&mut p, &mut max_t, &mut scale, &mut recovery, &mut o);
        parser.parse_args_or_exit();
    };
    p.resource_recovery_per_season = recovery * p.season_length_in_years;
    println!("# Initialization ...");
    // FIXME: There is something messed up in lat/long -> image pixels. It works
    // currently, but the names point out that I misunderstood something, eg.
    // with the layout of the tiff in memory or similar.

    let s = initialization(p, scale).unwrap();
    println!("Initialized");

    run(s, max_t, &o);
    Ok(())
}
