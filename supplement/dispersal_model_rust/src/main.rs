use model::argparse::parse_args;
use model::movementgraph::MovementGraph;
use model::{initialization, observation, run, Seasons};
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

    let mut p = model::submodels::parameters::Parameters {
        dispersal_graph,
        ..model::submodels::parameters::Parameters::default()
    };
    let mut max_t: Seasons = 20000;
    let mut scale = 1.0;
    let mut recovery = p.resource_recovery_per_season / p.season_length_in_years;

    let mut o = observation::Settings {
        log_every: 30, // Every 5 years, by default parameters
        log_patch_resources: 30,
        store_every: 600,
        statefile: "state.json".to_string(),
    };

    {
        let parser = parse_args(&mut p, &mut max_t, &mut scale, &mut recovery, &mut o);
        parser.parse_args_or_exit();
    };
    p.resource_recovery_per_season = recovery * p.season_length_in_years;
    println!("# Initialization ...");

    println!(" {:?}", o);
    println!(" {:?}", p);
    let s = initialization(p, scale);
    println!("Initialized");

    run(s, max_t, &o);
    Ok(())
}
