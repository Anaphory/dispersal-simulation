use crate::movementgraph::MovementGraph;
use crate::OneYearResources;
use serde_derive::{Deserialize, Serialize};

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct Parameters {
    pub resource_recovery_per_season: f64,
    pub culture_mutation_rate: f64,
    pub culture_dimensionality: usize,
    pub cooperation_threshold: u32,
    pub maximum_resources_one_adult_can_harvest: OneYearResources,
    pub evidence_needed: f64,
    pub payoff_std: f64,
    pub minimum_adaptation: f64,
    pub fight_deadliness: u32,
    pub enemy_discount: f64,

    pub season_length_in_years: f64,
    pub dispersal_graph: MovementGraph,
}

impl Default for Parameters {
    fn default() -> Parameters {
        Parameters {
            resource_recovery_per_season: 0.10,
            culture_mutation_rate: 6e-3,
            culture_dimensionality: 20,
            cooperation_threshold: 6,
            // In a previous iteration, there was no maximum of
            // resources accessible to one adult in a year. To be
            // consistent with that, we set the default of this new
            // parameter to a ludicrously high amount.
            maximum_resources_one_adult_can_harvest: OneYearResources::from(1e50),
            evidence_needed: 0.1,
            payoff_std: 0.1,
            minimum_adaptation: 0.5,
            fight_deadliness: (1 << 31),
            enemy_discount: (0.5_f64).powf(0.2),

            // I found that Kelly (2013), the guy who collected data
            // like Binford but maybe more faithfully to science and the
            // sources, has data on hunter-gatherer migrations, so I
            // could make a case for 6 migration steps/seasons per year,
            // and for estimating a maximum distance of each of these
            // steps. Together with a discussion I had with Peter and
            // Nico about two weeks ago, I managed to put something new
            // together. The downside is that with 6 seasons per year,
            // it also takes about a factor of 3 longer to simulate the
            // same number of years, so I'm only 250 years into the
            // simulation with this.
            season_length_in_years: 1. / 6.,
            dispersal_graph: MovementGraph::default(),
        }
    }
}
