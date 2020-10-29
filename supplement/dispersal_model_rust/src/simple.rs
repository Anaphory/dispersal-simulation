use model::hexgrid;
use model::movementgraph::MovementGraph;
use model::submodels::parameters::Parameters;
use model::argparse::parse_args;
use model::{run,State,KCal,initialization,HalfYears};
use std::collections::HashMap;
use std::fs;

fn main() -> Result<(), String> {
    let mut dispersal_graph = MovementGraph::default();

    let mut from = dispersal_graph.add_node(
        (0 as hexgrid::Index, 0., 0., vec![(0, 1.)].drain(..).collect()));
    for n in 1..23 {
        let to = dispersal_graph.add_node(
            (0 as hexgrid::Index, 0., -(n as f64), vec![(0, (2.0_f32).powi(n))].drain(..).collect()));
        dispersal_graph.add_edge(from, to, 1.);
        from = to
    }

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
        minimum_adaptation: 0.5,

        dispersal_graph: dispersal_graph,
    };
    let mut max_t: HalfYears = 2000;
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
