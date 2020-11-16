use crate::{Seasons, Parameters};
use crate::ecology::OneYearResources;

pub fn parse_args(p: &mut Parameters, max_t: &mut Seasons, resource_scale: &mut OneYearResources) {
    let mut parser = argparse::ArgumentParser::new();
    parser.set_description("Run a dispersal simulation");
    // TODO: Attach an ArgumentParser to the parameters object, so parameters can be set from the CLI.
    parser.refer(&mut p.attention_probability).add_option(
        &["--attention-probability"],
        argparse::Store,
        "attention probability",
    );
    parser.refer(&mut p.minimum_adaptation).add_option(
        &["--minimum-adaptation"],
        argparse::Store,
        "minimum adaptation of families to unknown ecoregions. 0: completely incapable, 1: fully adapted to all ecoregions.",
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
    parser.refer(resource_scale).add_option(
        &["--resources-scale"],
        argparse::Store,
        "Scale resources by this factor",
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
