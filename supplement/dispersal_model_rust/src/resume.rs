use ::argparse;
use model::*;
use std::fs;

fn main() -> Result<(), String> {
    let mut o = observation::ObservationSettings {
        log_every: 30, // Every 5 years, by default parameters
        log_patch_resources: 30,
        store_every: 600,
        statefile: "resumed.json".to_string(),
    };

    let mut statefile = "state.json".to_string();
    let mut end = 60000;
    {
        let mut parser = argparse::ArgumentParser::new();
        parser.set_description("Resume a dispersal simulation");
        parser.refer(&mut o.statefile).add_option(
            &["--statefile"],
            argparse::Store,
            "File to write state to",
        );
        parser.refer(&mut statefile).add_option(
            &["--resume-from"],
            argparse::Store,
            "File to read initial state from",
        );
        parser.refer(&mut o.log_every).add_option(
            &["--log-every"],
            argparse::Store,
            "period of logging, in seasons",
        );
        parser.refer(&mut end).add_option(
            &["--steps"],
            argparse::Store,
            "number of seasons to finish the simulation at (includes stored states)",
        );
        parser.refer(&mut o.store_every).add_option(
            &["--store-every"],
            argparse::Store,
            "period of state saving, in seasons",
        );
        parser.parse_args_or_exit();
    }

    let contents = match fs::read(statefile) {
        Ok(c) => c,
        Err(e) => return Err(e.to_string()),
    };

    let state: State = match bincode::deserialize(&contents) {
        Ok(c) => c,
        Err(e) => return Err(e.to_string()),
    };

    run(state, end, &o);
    Ok(())
}
