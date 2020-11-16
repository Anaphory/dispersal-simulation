use model::hexgrid;
use model::movementgraph::MovementGraph;
use model::submodels::parameters::Parameters;
use model::argparse::parse_args;
use model::ecology::OneYearResources;
use model::{run,initialization,Seasons};

fn main() -> Result<(), String> {
    let mut dispersal_graph = MovementGraph::default();

    let mut from = dispersal_graph.add_node(
        (0 as hexgrid::Index, 0., 0., vec![(0, 1.)].drain(..).collect()));
    for n in 1..23 {
        let to = dispersal_graph.add_node(
            (0 as hexgrid::Index, 0., -(n as f64), vec![(0, (2.0_f64).powi(n))].drain(..).collect()));
        dispersal_graph.add_edge(from, to, 1.);
        from = to
    }

    let mut p = Parameters::default();
    p.dispersal_graph = dispersal_graph;
    let mut max_t: Seasons = 2000;
    let mut scale = OneYearResources::from(1.0);
    parse_args(&mut p, &mut max_t, &mut scale);

    println!("# Initialization ...");
    let s = initialization(p, &scale).unwrap();
    println!("Initialized");

    let o = model::observation::ObservationSettings {
        log_every: 1,
    };
    run(s, max_t, &o);
    Ok(())
}
