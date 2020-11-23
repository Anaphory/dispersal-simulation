use model::hexgrid;
use model::movementgraph::MovementGraph;
use model::submodels::parameters::Parameters;
use model::argparse::parse_args;
use model::{run,initialization,Seasons};

fn main() -> Result<(), String> {
    let mut p = Parameters::default();
    let mut max_t: Seasons = 2000;
    let mut scale = 1.0;
    let mut spots = Vec::new();
    let mut recovery = p.resource_recovery_per_season / p.season_length_in_years;

    let mut o = model::observation::ObservationSettings {
        log_every: 1,
        log_patch_resources: 1,
    };

    {
        let mut parser = parse_args(&mut p, &mut max_t, &mut scale, &mut recovery, &mut o.log_every);
        parser.refer(&mut spots).add_option(
            &["--spot"],
            argparse::Collect,
            "create a spot with 2^SPOT population size",
        );
        parser.parse_args_or_exit();
    };
    p.resource_recovery_per_season = recovery * p.season_length_in_years;

    let mut dispersal_graph = MovementGraph::default();
    let mut from = None;
    for (i, n) in spots.iter().enumerate() {
        let to = dispersal_graph.add_node(
            (0 as hexgrid::Index, 0., -(i as f64), vec![(0, (2.0_f64).powi(*n))].drain(..).collect()));
        match from {
            None => (),
            Some(f) => {dispersal_graph.add_edge(f, to, 1.);}
        }
        from = Some(to)
    }
    p.dispersal_graph = dispersal_graph;

    println!("# Initialization ...");
    let s = initialization(p, scale).unwrap();
    println!("Initialized");

    run(s, max_t, &o);
    Ok(())
}
