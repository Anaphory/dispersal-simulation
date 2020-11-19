use crate::{Parameters, Seasons};

pub fn parse_args<'a>(
    p: &'a mut Parameters,
    max_t: &'a mut Seasons,
    resource_scale: &'a mut f64,
) -> argparse::ArgumentParser<'a> {
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
    parser.refer(&mut p.resource_density).add_option(
        &["--resource-scarcity"],
        argparse::Store,
        "Density of resources in a unit patch",
    );
    parser.refer(&mut p.storage_loss).add_option(
        &["--storage-loss"],
        argparse::Store,
        "Proportion of per-family stored resources getting lost per season",
    );
    parser.refer(&mut p.resource_recovery).add_option(
        &["--resource-recovery"],
        argparse::Store,
        "Proportion of resources to recover every season",
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
    parser
}
