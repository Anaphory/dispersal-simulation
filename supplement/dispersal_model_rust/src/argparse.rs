use crate::observation::ObservationSettings;
use crate::{Parameters, Seasons};
use std::str::FromStr;
use argparse::action::TypedAction;
use std::rc::Rc;
use std::cell::RefCell;
use argparse::action::{Action,IArgAction};
use argparse::action::ParseResult;
use argparse::action::ParseResult::{Parsed, Error};
use argparse::action::Action::Single;

pub struct StorePAsU32Action<'a> {
    pub cell: Rc<RefCell<&'a mut u32>>,
}

impl<'a> IArgAction for StorePAsU32Action<'a> {
    fn parse_arg(&self, arg: &str) -> ParseResult {
        match FromStr::from_str(arg) {
            Ok::<f64, _>(x) => {
                **self.cell.borrow_mut() = (x * ((1<<32) as f64)) as u32;
                return Parsed;
            }
            Err(_) => {
                return Error(format!("Bad value {}", arg));
            }
        }
    }
}

struct StorePAsU32;

impl TypedAction<u32> for StorePAsU32 {
    fn bind<'x>(&self, cell: Rc<RefCell<&'x mut u32>>) -> Action<'x> {
        return Single(Box::new(StorePAsU32Action { cell: cell }));
    }
}

pub fn parse_args<'a>(
    p: &'a mut Parameters,
    max_t: &'a mut Seasons,
    resource_scale: &'a mut f64,
    resource_recovery_per_year: &'a mut f64,
    observation: &'a mut ObservationSettings,
) -> argparse::ArgumentParser<'a> {
    let mut parser = argparse::ArgumentParser::new();
    parser.set_description("Run a dispersal simulation");
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
    parser.refer(&mut p.season_length_in_years).add_option(
        &["--season-length"],
        argparse::Store,
        "Length of a season/time step, in years",
    );
    parser
        .refer(&mut p.maximum_resources_one_adult_can_harvest)
        .add_option(
            &["--harvest-per-person-per-year"],
            argparse::Store,
            "Maximum resources one adult can harvest, per *year*",
        );
    parser.refer(&mut p.enemy_discount).add_option(
        &["--enemy-discount"],
        argparse::Store,
        "How much is the value of a location discounted for a single enemy (which is not countered by a friendly)? The default is to discount a location by 1/2 for every 5 enemies (so 0.87055056)",
    );
    parser.refer(&mut p.fight_deadliness).add_option(
        &["--fight-deadliness"],
        StorePAsU32,
        "Probability to die",
    );
    parser.refer(resource_recovery_per_year).add_option(
        &["--resource-recovery"],
        argparse::Store,
        "Proportion of resources to recover every year",
    );
    parser.refer(resource_scale).add_option(
        &["--resources-scale"],
        argparse::Store,
        "Scale resources by this factor",
    );
    parser.refer(&mut observation.store_every).add_option(
        &["--store-every"],
        argparse::Store,
        "period of state saving, in seasons",
    );
    parser.refer(&mut observation.log_every).add_option(
        &["--log-every"],
        argparse::Store,
        "period of logging, in seasons",
    );
    parser.refer(&mut observation.statefile).add_option(
        &["--statefile"],
        argparse::Store,
        "File to read state from and store back to",
    );
    parser.refer(max_t).add_option(
        &["--steps"],
        argparse::Store,
        "number of seasons to simulate",
    );

    println!("Starting…");
    parser
}
