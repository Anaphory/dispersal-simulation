use model::hexgrid;
use model::movementgraph::MovementGraph;
use model::submodels::parameters::Parameters;
use model::*;
use std::collections::HashMap;
use std::fs;

fn main() -> Result<(), String> {
    let contents = match fs::read("graph.bincode") {
        Ok(c) => c,
        Err(e) => return Err(e.to_string())
    };

    let dispersal_graph: MovementGraph = match bincode::deserialize(&contents) {
        Ok(c) => c,
        Err(e) => return Err(e.to_string())
    };

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

        dispersal_graph: dispersal_graph,
    };
    let mut max_t: HalfYears = 20000;
    parse_args(
        &mut p,
        &mut max_t,
    );
    println!("# Initialization ...");
    // FIXME: There is something messed up in lat/long -> image pixels. It works
    // currently, but the names point out that I misunderstood something, eg.
    // with the layout of the tiff in memory or similar.

    // p.dispersal_graph = ...
    let s: State = initialization(&p).unwrap();
    println!("Initialized");

    run(s, p, max_t);
    Ok(())
}

fn parse_args(
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
    println!("Starting…");
    parser.parse_args_or_exit();
}
