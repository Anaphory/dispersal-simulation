/*!
Model Description
=================

This model description follows the ODD (Overview, Design concept, Details)
protocol (Grimm et al., 2006; Grimm et al., 2010). The model description follows
the idea of literate programming (Knuth 1992) to the extent useful in Rust
source code – the actual model is generated from this file that documents the
model, but the source code is largely commented using natural-language
descriptions, not generated from them (as would be the case in a literal
program).

# 1. Purpose

The dispersal model generates a deep phylogeny of hunter-gatherer cultures based
on culture-mediated cooperation and resource-driven migration. It is a
demographic migration model in which the areal distribution of languages is an
emergent property, not an imposed structure. The model is designed with
extension to more concrete research questions in mind. In the current, first
stage, the purpose of the model is to investigate how languages disperse and
split, driven only by the necessary interactions between humans. The particular
research question is whether distinct cultural areas, with reasonably sharp
borders, can emerge from only a simple feedback loop between conditional
cooperation, cultural assimilation, and genetic drift.

If that is the case, the summary statistics of a phylogeny (in particular
diversification rates) can later be compared to values known from language
evolution. The model is structured to be easily applied to study the history of
the settlement of the Americas at a later time. It would require paleoclimate
data and interpretations of past ecoregions to produce results that can be
compared to that history.

 */

// TODO: Dear Rustacean, I know that my use of documentation comments is
// hazardous, because they land in the generated documentation in a different
// order, or sometimes attached to things that are only accidentally coming afterwards.

// Load useful modules

use serde_derive::{Deserialize, Serialize};
use dashmap::DashMap;

use std::fs::File;
use std::io::prelude::*;
use rand::prelude::*;
use std::collections::HashMap;

use rayon::prelude::*;

pub mod argparse;
pub mod ecology;
mod debug;

use submodels::parameters::Parameters;

/**

# 2. Entities, state variables, and scales

The model consists of agents interacting on a network of habitable patches in
discrete time. One time step is supposed to model a season with a duration of
half a year. Some of the underlying data is using seconds as base unit and needs
to be translated.

 */
pub type Seasons = u32;
const SECONDS_PER_YEAR: f64 = 365.24219 * 24. * 60. * 60.;

/**
Whereever possible, resources are measured in terms of the resources one adult
consumes in a year, abbreviated OYR.
 */
use ecology::OneYearResources;

/**
## 2.1 Grid and Patches

The geography of the simulation is described by a directed weighted graph. Each
node is a hexagon from a hexagonal equal-area discrete global grid, each edge
has a weight representing the travel time (by foot, or by simple boat along
rivers and coasts) in seconds. Each node has an arbitrary unique numerical
identifier. This movement graph is constructed before the initialization from
pre-processed real geographical data.

 */
pub mod movementgraph;
pub type NodeId = petgraph::graph::NodeIndex<usize>;

/**
Each node contains one or more ecoregions, according to the geographical
distribution of ecoregions in the real world. For each ecoregion, the patch has
a maximum and current availability of resources, measured in OYR available over
the course of one season. These numbers are estimated from population densities,
adapting the estimation of population density per biome from [@tallavaara] to
obtain a raster of population densities, which together with the shapes of the
hexagons gives population counts per ecoregion area in km².

In the current state of the model, the resource maximum is constant throughout
the simulation, but future work might add seasonality (with every other time
step corresponding to the more plentiful half of the year and the others to the
scarcer half) or random variation. By nature, the simulation invites the
extension to include paleoclimate data, such as recently published in
[@beyer2020highresolution] but that only becomes relevant once it shows
fundamentally reasonable dynamics.

https://www.nature.com/articles/s41597-020-0552-1

 */
#[derive(Serialize, Deserialize, Clone)]
pub struct Patch {
    /// For every local ecoregion, the tuple of resources currently accessible
    /// in this patch, and the maximum possible resources available in a time
    /// step, both in OYR.
    resources: HashMap<usize, (OneYearResources, OneYearResources)>,
}

/**
## 2.2 Families

The main decision-making agents of the simulation are families. Families can
migrate between nodes and form links to other families in the context of
cooperation to extract resources.

 */

#[derive(Serialize, Deserialize)]
pub struct Family {
    /// The agent's history of decendence, also serving as unique ID.
    pub descendence: String,
    /// The effective size of the family in number of adults. One adult is
    /// assumed to consume the same amount of food/energy and to contribute the
    /// same labor to foraging as any other adult.
    pub effective_size: usize,
    /// The family's shared culture.
    pub culture: Culture,
    /// The current location of the agent.
    pub location: NodeId,
    /// The amount of stored resources, in OYR, the family has access to
    /// without going foraging.
    pub stored_resources: OneYearResources,
    /// Adaptation to local ecoregions is represented as a vector with values
    /// between 0.5 (completely unknown) and 1.0 (perfectly familiar) for any
    /// ecoregion the family might encounter.
    pub adaptation: ecology::Ecovector,
    /// The number of seasons to wait until the next child (reset to 2 when
    /// starvation happens)
    pub seasons_till_next_child: Seasons,

    /// For bookkeeping purposes (eg. generating descendant's ‘descendence’
    /// values), keep track of the number of offspring families this family has
    /// spawned so far.
    pub number_offspring: u16,
    /// The previous locations of the agent. This is useful for some bits of
    /// analysis.
    pub history: Vec<NodeId>,
    /// A bookkeeping quantity. Instead of mutation happening in each time step
    /// with a tiny probability, the same distribution of times between
    /// mutations is generated by drawing a number of season from a geometric
    /// distribution and counting it down by one each time step, which is useful
    /// because random number generation is computationally expensive.
    pub seasons_till_next_mutation: Option<Seasons>,
}

impl PartialEq for Family {
    fn eq(&self, other: &Self) -> bool {
        self.descendence == other.descendence
    }
}

/**

Families in the same location with compatible cultures can band together to
share resources between different families, and to compete against other bands.
The cooperative bands formed by cooperating families are higher-level agents
created ad-hoc in each time step. They do not persist or have effect beyond a
single time step. Resource exploitation happens at the level of the cooperative
group and is distributed to the individual families after the fact.

 */
use collectives::Cooperative;

/**
## 2.3 Cultures

Every family has a culture. These are very abstract and vastly simplified, due
to the lack of quantitative data on cultural evolution in a framework comparable
to the one used for this study. Based on the need to have a data structure that
supports random drift, a low but positive chance of back-mutation, and a large
number of equally (dis-)similar cultures, we describe culture using a binary
vector. Computers make natural use of the binary representation integer numbers,
so a binary vector is equivalent to an unsigned integer of sufficient size,
which is faster to use in computations and more efficient to store.

 */
#[derive(Ord, PartialOrd, Eq, PartialEq, Hash, Clone, Copy, Serialize, Deserialize)]
pub struct Culture {
    binary_representation: u64,
}

impl std::fmt::Binary for Culture {
    fn fmt<'a>(&self, f: &mut std::fmt::Formatter<'a>) -> Result<(), std::fmt::Error> {
        self.binary_representation.fmt(f)
    }
}

impl From<u64> for Culture {
    fn from(v: u64) -> Culture {
        Culture {
            binary_representation: v,
        }
    }
}

/**
## 2.4 State

The associations between nodes and patches and between nodes and the families
located there (stored in the Family's `location`) are the core of the model
state. The state also tracks the time, measured in time steps corresponding to
one season each, since the start of the simulation, and stores a copy of
the model parameters.

 */
#[derive(Debug, Serialize, Deserialize)]
pub struct State {
    /// The patches of the model, indexed by node ID in a graph. This set is
    /// fixed, although the properties of the patches may change with time.
    patches: DashMap<NodeId, Patch>,
    /// The agents (families) currently active in the simulation. This vector
    /// changes over time as new families are added and old families die out.
    families: Vec<Family>,
    /// The current time step, in half years since start of the simulation.
    pub t: Seasons,
    /// The parameters of the model
    p: Parameters,
}

/**
# 3. Process overview and scheduling

>>> Who (i.e., what entity) does what, and in what order? When are state
>>> variables updated? How is time modeled, as discrete steps or as a continuum
>>> over which both continuous processes and discrete events can occur? Except
>>> for very simple schedules, one should use pseudo-code to describe the
>>> schedule in every detail, so that the model can be re-implemented from this
>>> code. Ideally, the pseudo-code corresponds fully to the actual code used in
>>> the program implementing the ABM.

The model progresses in discrete time steps, each corresponding to half a year
of simulated time. The entire simulation consists of repeating the step to
simulate 15000 years.

The structure of a single time step consist of two parts, as follows.

 */
fn step(
    families: &mut Vec<Family>,
    patches: &mut DashMap<NodeId, Patch>,
    p: &Parameters,
    t: Seasons,
    o: &observation::ObservationSettings,
) {
    step_part_1(families, patches, p);
    step_part_2(families, patches, p, o, t);
}
/**
The first part focuses on the individual families, which shrink, grow, die,
split, and move. It constructs the mapping of families at the end of the season,
grouped by their location after potential moves. Because the movement of a
family depends on the distribution of the families at he start of the season and
not at the time of movement, it can happen entirely in parallel.

 */
fn step_part_1(families: &mut Vec<Family>, patches: &mut DashMap<NodeId, Patch>, p: &Parameters) {
    let cc: DashMap<NodeId, usize> = DashMap::new();
    for f in families.iter() {
        *cc.entry(f.location).or_insert(0) += f.effective_size;
    }

    families.par_iter_mut().for_each(|mut family| {
        if submodels::family_lifecycle::use_resources_and_maybe_shrink(
            &mut family.effective_size,
            &mut family.stored_resources,
            p,
        ) {
            family.seasons_till_next_child = std::cmp::max(family.seasons_till_next_child, 2)
        }

        submodels::family_lifecycle::maybe_grow(&mut family, p.season_length_in_years);
    });

    families.retain(|family| submodels::family_lifecycle::can_survive(&family));

    let mut children =
        submodels::family_lifecycle::procreate_and_migrate(families, p, patches, &cc);
    families.extend(children.drain(..));
}

/**
The second part focusses on the patches. The resources of a patch are updated
according to the families exploiting them over the season. This is described in
detail in Submodule 7.5. Everything here happens locally to a patch, with no
external interaction, so this can be done in parallel. After exploitation,
patches recover advance to the next season according to Submodule 7.6. This
concludes a time step.

 */
fn step_part_2(
    families: &mut Vec<Family>,
    patches: &mut DashMap<NodeId, Patch>,
    p: &Parameters,
    o: &observation::ObservationSettings,
    t: Seasons,
) {
    let families_by_location: DashMap<NodeId, Vec<&mut Family>> = DashMap::new();
    families.par_iter_mut().for_each(|f| {
        families_by_location
            .entry(f.location)
            .or_insert_with(Vec::default)
            .push(f)
    });

    let cultures_by_location: HashMap<NodeId, HashMap<Culture, usize>> = families_by_location
        // TODO: Change when DashMap directly supports par_iter
        .into_iter()
        .par_bridge()
        .map(|(patch_id, families)| {
            let mut cc = HashMap::new();

            let mut groups = crate::collectives::cooperatives(families, p);

            for group in groups.iter_mut() {
                let group_culture = submodels::ecology::adjust_culture(group, p);
                cc.insert(
                    group_culture,
                    group.families.iter().map(|f| f.effective_size).sum(),
                );
            }

            match patches.get_mut(&patch_id) {
                None => {}
                Some(patch) => {
                    if (o.log_patch_resources > 0) && (t % o.log_patch_resources == 0) {
                        // println!("Patch at {:?}: {:?}", patch_id, patch.resources)
                    }
                    submodels::ecology::exploit_patch(groups, patch, p);
                }
            };

            (patch_id, cc)
        })
        .collect();

    if (o.log_every > 0) && (t % o.log_every == 0) {
        println!("t: {:}", t);
        observation::print_population_by_location(&cultures_by_location, p);
        println!("F: {:?}", families.get(0));
    }
    if (o.log_gdcd > 0) && (t % o.log_gdcd == 0) {
        println!("t: {:}", t);
        observation::print_gd_cd(&cultures_by_location, p);
    }
}

/**
# 4. Design concepts

Under the ODD protocol, the design principles largely fall into questions. Where
my answer indicates an invariant property of the simulation, I provide a test
function that checks that invariance if possible.

 */
mod concepts {}

/**
## 4.1 Basic priciples

> Which general concepts, theories, hypotheses, or modeling approaches are
> underlying the model’s design?

According to its purpose of providing a model for language dispersal and split,
our model draws on existing publications looking at for cultures changing in
time and space, focussing on taking a bottom-up approach on languages splitting.
A major useful reference is the PSMED model (delcastillo2013modeling,
barcelo2013psmed, barcelo2015simulating), which includes two drivers of cultural
change in the agents: Drift and assimilation. In isolation, agents' cultures
undergo drift. Where social interactions happen, drift is counterbalanced by
cultural assimilation for agents that cooperate with each other. In the case of
that model, agents interact with other agents in a limited range and based on
need. We posit a similar basic relationship between *cooperation driving
cultural assimilation, and cultural similarity conditioning cooperation*. Due to
the larger scope of our model, however, we avoid tracking networks of pairwise
interactions and instead assume that cooperation happens in patches,
unconditionally between all individuals of compatible culture in that patch.

This interaction between cooperation and culture is the main addition to an
otherwise *demographic dispersal model*. The migratory pattern we will see are
not imposed by any geographical structure, buy by the carrying capacities of the
individual patches. Difficult terrain (mountains, seas) is, to the model, just
terrain with a low carrying capacity, but does not a priori prevent agents from
crossing it if they can reach a patch on the other side that is habitable. Due
to a limited interaction range, however, larger uninhabitable areas do provide
an impediment to movement, because agents naturally avoid moving into them if
better options exist.

A major necessary feature of the dispersal component for the culture component
is to allow *group fission and fusion*, to allow culture splits to emerge from
low-level agent interactions, but to also drive incentives for agents to
congregate, cooperate, and assimilate their cultures. As such, the model draws
heavily on the agent-based model by crema2014simulation, crema2015modeling, with
minor variants.

FIXME: Given that statement, I should really construct an implementation of
Crema's model compatible with my setup, where the differences are explicitly
visible as model parameters that I use differently from them. My patch topology
is different, what else?

 */

/**
One of the deviations from Crema's model, and also from PSMED, is the geography
underlying the simulation. Crema and Barceló use a quadratic grid with arbitrary
fixed or randomly generated resources. A long-term goal for our model it
applicability to real-world language dispersal processes. A hexagonal grid has
some advantages over a quadratic grid on the plane, in that all hexagons within
a fixed number of steps from a focus are a much better approximation of a circle
than all squares within a fixed number of steps from a central square. Hexagons
also have nicer properties for delineating boundaries, which need some
approximation for square grids (cf. kuijper2004detecting, which is not directly
on finding boundaries, but solving that problem nonetheless). On a sphere such
as Earth, using a hexagonal grid has the additional advantage of reducing
distortions. Using hexagonal grids in continental or global dispersal
simulations is well-established (gavin2017processbased, callegari2013agentbased).

Our model is thus an agent-based model of the emergence of cultural areas from
fission-fusion dynamics, where local cooperation is the main driver to prevent
cultures from unconstrained evolutionary drift.

> How were they taken into account?

> Are they used at the level of submodels, or is their scope the system level?

> Will the model provide insights about the basic principles themselves, i.e.,
> their scope, their usefulness in real-world scenarios, validation, or
> modification?

> Does the model use new, or previously developed, theory for agent traits from
> which system dynamics emerge?

 */

mod basic_principles {}

/**
## 4.2 Emergence

> What key results or outputs of the model are modeled as emerging from the
> adaptive traits, or behaviors, of individuals? In other words, what model
> results are expected to vary in complex and perhaps unpredictable ways when
> particular characteristics of individuals or their environment change?

*/
mod emergence {
    use crate::*;

    /**

    The main emergent property will be culturally somewhat uniform territories. We
        expect that the interplay of migration, cooperation and cultural similarity
        leads to regions sharing similar cultures, with noticeable boundaries. This
        means that the plot of cultural distances vs. geographical distances should
        therefore show small cultural distances for small geographical distances. There
        should be a critical geographical distance of cohesion where the distribution of
        cultural distances becomes bimodal, with one mode being lower than the
        cooperation threshold and one mode above the cooperation threshold. For large
        geographical distances, the cultural distances should be high, but with a big
        variance. The critical geographical distance is likely to depend on the region,
        being larger in more marginal environments where migration is more frequent.

     */
    pub fn cultural_distance_by_geographical_distance(
        cultures_by_location: &HashMap<NodeId, HashMap<Culture, usize>>,
        max_geo_distance: i32,
        p: &Parameters,
    ) -> HashMap<i32, HashMap<u32, usize>> {
        let mut gd_cd: HashMap<i32, HashMap<u32, usize>> = HashMap::new();
        for (l1, cs1) in cultures_by_location {
            let l1h3 = p.dispersal_graph.node_weight(*l1).unwrap().0;
            for (l2, cs2) in cultures_by_location {
                let l2h3 = p.dispersal_graph[*l2].0;
                let gd = hexgrid::hex_distance(l1h3, l2h3);
                if gd > max_geo_distance {
                    continue;
                }
                for (c1, count1) in cs1 {
                    for (c2, count2) in cs2 {
                        if l1 == l2 && c1 == c2 {
                            *gd_cd
                                .entry(gd)
                                .or_insert_with(HashMap::new)
                                .entry(0)
                                .or_insert(0) += count1 * (count2 - 1) / 2;
                            break;
                        }
                        let cd = submodels::culture::distance(*c1, *c2);
                        *gd_cd
                            .entry(gd)
                            .or_insert_with(HashMap::new)
                            .entry(cd)
                            .or_insert(0) += count1 * count2;
                    }
                }
                if l1 == l2 {
                    break;
                }
            }
        }
        gd_cd
    }

    /**

    Underlying the model is a fission-fusiom model of group dynamics
        (crema2014simulation), so we expect a similar analysis to apply to our
        model. Crema observes actual fission-fusion cycles only for a small part of
        the parameter space (“Primate systems emerge only temporarily, either as
        part of a limit-cycle equilibrium, or as a short-term transition from a
        convex equilibrium. In either case, they require some level of system
        connectivity, defined here by the spatial range of interaction (h), the
        frequency of decision-making (z), and the sample proportion of observed
        neighbour agents (k)”), so it is not clear whether we should expect them for
        our model.

    > Are there other results that are more tightly imposed by model rules and
    > hence less dependent on what individuals do, and hence ‘built in’ rather
    > than emergent results?

    Cooperation in this model is an entirely imposed feature. It is an important
        question in theoretical biology and related fields under what circumstances
        cooperation can prevail over defectors and other free-riders, and it may
        even affect the dispersal of languages (though more likely on the level of
        collectives, where the interactions between groups of different cultures
        range from assimilation and language shift all the way to war and lethal
        violence). But overall, human societies seem to be quite good at maintaining
        large-scale cooperation, so we consider the question irrelevant for the
        purpsoses of the present model. Cooperation is thus not a decision of the
        agents, but entirely determined by their cultures.

     */

    pub fn similar_culture(c1: Culture, c2: Culture, cooperation_threshold: u32) -> bool {
        submodels::culture::distance(c1, c2) < cooperation_threshold
    }
}

/**
## 4.3 Adaptation

> What adaptive traits do the individuals have? What rules do they have for
> making decisions or changing behavior in response to changes in themselves or
> their environment? Do these traits explicitly seek to increase some measure of
> individual success regarding its objectives? Or do they instead simply cause
> individuals to reproduce observed behaviors that are implicitly assumed to
> indirectly convey success or fitness?

The only trait agents have control over is their location. Agents optimize
(within the limits of their knowledge) their location to increase the amount of
resources they gather. The resources gathered indirectly contribute to future
success: A certain minimum of gathered resources is necessary to procreate, but
extremely large amounts gathered soon do not increase the immediate rate of
procreation and can in extreme cases (overexploitation of local resources) even
be detrimental in the long run.
 */

/**
Agents also adapt to local environment: TODO

 */
mod adaptation {
    use crate::*;

    pub fn decide_on_moving<'a, KD>(
        family: &'a Family,
        known_destinations: KD,
        patches: &'a DashMap<NodeId, Patch>,
        avoid_stay: bool,
        p: &Parameters,
        population: &DashMap<NodeId, usize>,
    ) -> (NodeId, OneYearResources)
    where
        KD: Iterator<Item = (NodeId, f64)>,
    {
        let evidence = OneYearResources::from(p.season_length_in_years) * p.evidence_needed;

        // Triple: opportunity cost from not foraging, and cost from doing likely
        // heavy labor instead.
        let cost_per_walking_second: OneYearResources =
            OneYearResources::from(1.) / SECONDS_PER_YEAR * 3.;

        let destination_expectation = known_destinations.filter_map(|(i, d)| {
            if family.location == i && avoid_stay {
                return None;
            }
            let movement_cost = cost_per_walking_second * d * family.effective_size as f64;
            let pop = match population.get(&i) {
                None => 0,
                Some(k) => *k,
            } + (if family.location == i {
                0
            } else {
                family.effective_size
            });
            let now = expected_quality(i, p, patches, pop);

            Some((i, (now, movement_cost)))
        });
        objectives::best_location(destination_expectation, evidence)
            .unwrap_or((family.location, OneYearResources::from(0.0)))
    }

    fn expected_quality(
        i: NodeId,
        p: &Parameters,
        patches: &DashMap<NodeId, Patch>,
        population: usize,
    ) -> OneYearResources {
        // Minimal knowledge about quality of a patch: Its current plus max resources.
        let mut rng = rand::thread_rng();
        let r = patches.get(&i).unwrap();
        let q: Vec<_> = r.resources.values().collect();
        let perhead = 1. / population as f64;
        let now = q
            .iter()
            .map(|(res, max_res)| *res + *max_res * p.resource_recovery_per_season)
            .sum::<OneYearResources>()
            * perhead
            * rng.gen_range(0., 1.);
        // println!("Option: {:?}, with current resources {:?} (per head: {:?}, for {:}) at {:?}", i, (now + movement_cost) / perhead, now, population, movement_cost);
        now
    }
}

/**
## 4.4 Objectives

> If adaptive traits explicitly act to increase some measure of the individual’s
> success at meeting some objective, what exactly is that objective and how is
> it measured?

> When individuals make decisions by ranking alternatives, what criteria do they
> use?

*/
mod objectives {
    use crate::*;

    /**
    Agents choose a best location to move to, where ‘best’ means a maximum
    expected resource gain. If multiple locations are equally good, take one of those at random.
    */
    // That is, this function omputes the argmax of
    // `expected_resources_from_patch`, drawing at random between equal options.
    pub fn best_location<Pair, I>(
        kd: Pair,
        precision: OneYearResources,
    ) -> Option<(I, OneYearResources)>
    where
        Pair: Iterator<Item = (I, (OneYearResources, OneYearResources))>,
        I: std::fmt::Debug,
    {
        let mut rng = rand::thread_rng();

        // This variable `n_best_before` is used to randomly draw between
        // several equally-optimal options. It counts the number of best options
        // encountered so far.
        let mut n_best_before = 0;
        let mut n_best_short_before = 0;
        let mut target: Option<(I, OneYearResources)> = None;
        let mut max_gain = OneYearResources::from(0.0);
        let mut max_short_term_gain = OneYearResources::from(0.0);

        for (location, (expected_shortterm_gain, travel_cost)) in kd {
            // println!("{:?}: short {:?} move {:?}", location, expected_shortterm_gain, travel_cost);
            if expected_shortterm_gain > max_gain + precision {
                target = Some((location, travel_cost));
                max_gain = expected_shortterm_gain;
                n_best_before = 0;
            } else if expected_shortterm_gain >= max_gain - precision {
                n_best_before += 1;
                if rng.gen_range(0, n_best_before + 1) < n_best_before {
                    continue;
                } else {
                    target = Some((location, travel_cost));
                }
            }
        }
        // println!("Moving to {:?}", target);
        target
    }
}
/**
## 4.5 Learning

> Many individuals or agents (but also organizations and institutions) change
> their adaptive traits over time as a consequence of their experience? If so,
> how?

There is a minor effect of experience on location choice in that a family knows
nearby locations only with a small random chance, but they always have knowledge
about the locations they were at in the last 4 years (8 time steps).

TODO: Put that 4 in a static variable and put a doctest here that confirms it.

 */

mod learning {}

/**
## 4.6 Prediction

> Prediction is fundamental to successful decision-making; if an agent’s
> adaptive traits or learning procedures are based on estimating future
> consequences of decisions, how do agents predict the future conditions (either
> environmental or internal) they will experience? If appropriate, what internal
> models are agents assumed to use to estimate future conditions or consequences
> of their decisions? What tacit or hidden predictions are implied in these
> internal model assumptions?

Agents predict and hypothesize their payoff, in terms of extracted resources, at
the end of the time step and decide accordingly. The actual estimation procedure
is modeled on a very abstract level: The payoff an agent predicts for themselves
in a location is internally computed using the same function that genererates
the actual payoff at the exploitation step of the simulation, without random
effects, thus giving the expectation instead of the actual value.

The prediction is based on the current state of the model, with very little
extrapolation whatsoever. An agent only adds their own effective size to the
current distribution of individuals in the target location. Agents do not
account for the following moves of other agents, and even less so for growth,
split, or other effects to themselves or others in future time steps.

 */
mod prediction {}

/**
## 4.7 Sensing

> What internal and environmental state variables are individuals assumed to
> sense and consider in their decisions? What state variables of which other
> individuals and entities can an individual perceive; for example, signals that
> another individual may intentionally or unintentionally send? Sensing is often
> assumed to be local, but can happen through networks or can even be assumed to
> be global (e.g., a forager on one site sensing the resource levels of all
> other sites it could move to). If agents sense each other through social
> networks, is the structure of the network imposed or emergent? Are the
> mechanisms by which agents obtain information modeled explicitly, or are
> individuals simply assumed to know these variables?

We take this a bit backwards: Kelly (2013), Table 4.1 (reproduced in the
supplementary material) lists mobility data for 70 hunter-gatherer groups,
including figueres for the yearly territory of a group where that number could
be inferred. Three quarters of the listed territories are about 2000 km² or
less. In our model, this corresponds to roughly 5 patches (each with an area of
ca. 400 km²) being visited per year. To give some leeway, we permit our agents
up to 6 moves per year, so a season corresponds to two months. Three quarters of
the total yearly distance the groups move are reported to be 320 km or less, so
in this simulation the maximum distance moved per season is about 53.3 km. The
simulation is centered around distances measured in time units (seconds,
seasons, years). In all but mangroves, the terrain speed factor is greater than
one (see `supplement/distances/ecoreginos.py`). If humans generally prefer 2%
downhill slopes, they have to encounter about the same amount in uphill slopes,
which are navigated at about 0.7740708353271509 m/s. At that speed, 53.3 km
correspond do about 68899 s or 19 hours. We therefore give individuals a maximum
“sensing” range of 68899 seconds. Individuals will know about patches within
this range from their current patch, but not about patches further away.

 */

mod sensing {
    use crate::*;

    /// Individuals know about nearby locations. The exploration is not
    /// explicitly modelled.

    pub fn nearby_locations(location: NodeId, p: &Parameters) -> HashMap<NodeId, f64> {
        movementgraph::bounded_dijkstra(
            &p.dispersal_graph,
            location,
            68899., // NORM * 12., //4 complete days, or 12 days of travel
            |e| *petgraph::graph::EdgeReference::weight(&e),
        )
    }
}

/**
## 4.8 Interaction

> What kinds of interactions among agents are assumed? Are there direct
> interactions in which individuals encounter and affect others, or are
> interactions indirect, e.g., via competition for a medi- ating resource? If
> the interactions involve communication, how are such communications
> represented?

Agents compete for limited resources. There is a maximum of resources that can
be extracted from a patch, and the extracted resources are distributed among all
agents in that patch. This distribution is not completely even: Members of a
bigger group of cooperators competing with a smaller cooperative get an
over-proportional part of the total extracted, measured by their effective size.

 */
mod interaction {
    /**

    A single human can forage a certain amount, but by the assumptions of the model,
        two cooperating foragers gain more resources together than they would
        individually. The effective labor returned by this function is the factor by
        which 2 cooperating foragers are better than 2 separate foragers.

    This formula follows Crema (2014), with the caveat that Crema computes payoffs,
        whereas we compute a multiplicative factor (effective labor) which is multiplied
        with the standard resource gain to give the acutal resource payoffs. Given that
        Crema's payoffs have a mean of μ=10, we adjust our effective labor by the
        opposite factor.

    TODO: I actually dislike this weakly motivated k → k + m k ^ (α + 1) quite a
        bit, so it would be nice to take a different formula, and to cross-check how
        Crema's results change for that different formula.

     */
    pub fn x() {}
}
/**

*/

/**
## 4.9 Stochasticity

> What processes are modeled by assuming they are random or partly random? Is
> stochasticity used, for example, to reproduce variability in processes for
> which it is unimportant to model the actual causes of the variability? Is it
> used to cause model events or behaviors to occur with a specified frequency?

TODO: Is there a way to get rust to aggregate all functions that use it? Can I
maybe somehow abuse the trait system to all those locations show up here, with
comments?

 */

/**
## 4.10 Collectives

> Do the individuals form or belong to aggregations that affect, and are
> affected by, the individuals? How are collectives represented? Is a particular
> collective an emergent property of the individuals, such as a flock of birds
> that assembles as a result of individual behaviors, or is the collective
> simply a definition by the modeler, such as the set of individuals with
> certain properties, defined as a separate kind of entity with its own state
> variables and traits?

 */
pub mod collectives {
    use crate::*;
    pub struct Cooperative<'a> {
        pub families: Vec<&'a mut Family>,
        pub culture: Culture,
    }

    pub fn cooperatives<'a>(
        mut families_in_this_location: Vec<&'a mut Family>,
        p: &Parameters,
    ) -> Vec<Cooperative<'a>> {
        let mut rng = rand::thread_rng();
        let mut groups: Vec<Cooperative<'a>> = vec![];

        for family in families_in_this_location.drain(..) {
            let mut joined_group: Option<&mut Cooperative<'a>> = None;
            for group in groups.iter_mut() {
                let mut join = true;
                for other_family in group.families.iter() {
                    if !emergence::similar_culture(
                        family.culture,
                        other_family.culture,
                        p.cooperation_threshold,
                    ) {
                        join = false;
                        break;
                    }
                }
                if join {
                    joined_group = Some(group);
                    break;
                }
            }
            let _size = family.effective_size as f64;
            match joined_group {
                None => {
                    let group = Cooperative {
                        culture: family.culture,
                        families: vec![family],
                    };
                    groups.push(group);
                }
                Some(group) => {
                    if rng.gen_range(0, group.families.len()) == group.families.len() {
                        group.culture = family.culture
                    }
                    group.families.push(family);
                }
            }
        }
        groups
    }
}

/*
# 4.11 Observation

> What data are collected from the ABM for testing, understanding, and analyzing
> it, and how and when are they collected? Are all output data freely used, or
> are only certain data sampled and used, to imitate what can be observed in an
> empirical study (‘virtual ecologist’)?

 */
pub mod observation {
    use crate::*;
    /**
    Do get an overview over the flow of human migration captured by the model,
    we collect the population count, i.e. the sum of all families' effective
    sizes, for each spot in each time step.
    */
    pub fn print_population_by_location(
        cultures_by_location: &HashMap<NodeId, HashMap<Culture, usize>>,
        p: &Parameters,
    ) {
        println!(
            "POPULATION: {:?}",
            cultures_by_location
                .par_iter()
                .map(|(k, v)| {
                    let node = &p.dispersal_graph[*k];
                    (node.1, node.2, v)
                })
                .collect::<Vec<(f64, f64, _)>>()
        );
    }

    pub struct ObservationSettings {
        pub log_every: Seasons,
        pub log_gdcd: Seasons,
        pub store_every: Seasons,
        pub log_patch_resources: Seasons,
        pub statefile: String,
    }

    /**
    The major focus of the model, however, is on the cultural areas and
    networks. To asses the main outcome of the model, we collect the geographic
    vs. cultural distances for every pair of individuals within 5 years of each other. Due to the
    computational effort and high autocorrelation of the outcomes, we collect
    this data only about once per generation, i.e. every 60 time steps.
    */
    pub fn print_gd_cd(
        cultures_by_location: &HashMap<NodeId, HashMap<Culture, usize>>,
        p: &Parameters,
    ) {
        let gd_cd: HashMap<i32, HashMap<u32, usize>> =
            emergence::cultural_distance_by_geographical_distance(
                cultures_by_location,
                5 * 5 * 2,
                p,
            );
        println!("GD_CD: {:?}", gd_cd)
    }
}

/**
5. Initialization

> What is the initial state of the model world, i.e., at time $t = 0$ of a
> simulation run?

> In detail, how many entities of what type are there initially, and what are
> the exact values of their state variables (or how were they set
> stochastically)?

> Is initialization always the same, or is it allowed to vary among simulations?

> Are the initial values chosen arbitrarily or based on data? References to
> those data should be provided.

*/
pub mod hexgrid;

fn very_coarse_dist(x1: f64, y1: f64, x2: f64, y2: f64) -> f64 {
    (x1 - x2).abs() + (y1 - y2).abs()
}

pub fn initialization(p: Parameters, scale: f64) -> Option<State> {
    let graph = &p.dispersal_graph;

    let mut m_start1: Option<NodeId> = None;
    let mut m_start2: Option<NodeId> = None;
    let mut start1d = f64::INFINITY;
    let mut start2d = f64::INFINITY;

    for i in graph.node_indices() {
        let longitude = graph[i].1;
        let latitude = graph[i].2;
        let d1 = very_coarse_dist(longitude, latitude, -159.873, 65.613);
        let d2 = very_coarse_dist(longitude, latitude, -158.2718, 60.8071);
        if d1 < start1d {
            m_start1 = Some(i);
            start1d = d1;
        }
        if d2 < start2d {
            m_start2 = Some(i);
            start2d = d2;
        }
    }
    let start1 = m_start1.unwrap();
    let start2 = m_start2.unwrap();

    println!("Starts: {:}, {:}", start1.index(), start2.index());

    let mut patches: HashMap<NodeId, Option<Patch>> = HashMap::new();
    let mut new_patches: std::collections::BinaryHeap<NodeId> = std::collections::BinaryHeap::new();
    new_patches.push(start1);
    new_patches.push(start2);
    loop {
        let next = match new_patches.pop() {
            None => break,
            Some(n) => n,
        };
        if patches.contains_key(&next) {
            continue;
        }
        let geo = &graph[next];
        let longitude = geo.1;
        let latitude = geo.2;
        let ecoregions = &geo.3;

        print!(
            "Resources at {:} ({:}, {:})",
            next.index(),
            longitude,
            latitude
        );

        patches.insert(
            next,
            Some(Patch {
                resources: ecoregions
                    .iter()
                    .map(|(&ecoregion, &popcap)| {
                        // This over-estimates the carrying capacity of
                        // transition regions: For two ecoregions, each with
                        // carrying capacity 16, and cooperation gain 1, this
                        // gives twice 4 resources, which when harvested give 8
                        // resources together, which under cooperation will lead
                        // to a capacity of 8²=64>32=16+16. &>Transitory regions
                        // should be small, important for the movement between
                        // the ecoregions, and are otherwise disadvantaged (see
                        // the adaptation formula), so this is hopefully not too
                        // bad.
                        let q = popcap * scale / p.resource_recovery_per_season;
                        let res = OneYearResources::from(q);
                        let res_max = OneYearResources::from(q);
                        print!(
                            ", Filling ecoregion {:} with popcap {:} with {:?}",
                            ecoregion, popcap, res
                        );
                        (ecoregion, (res, res_max))
                    })
                    .collect(),
            }),
        );
        println!(".");
        for q in graph.neighbors(next) {
            if patches.contains_key(&q) {
                continue;
            }
            new_patches.push(q);
        }
    }

    Some(State {
        patches: patches
            .drain()
            .filter_map(|(i, p)| match p {
                None => None,
                Some(q) => Some((i, q)),
            })
            .collect(),
        families: vec![
            Family {
                descendence: String::from("A"),
                location: start1,
                history: vec![],
                seasons_till_next_child: 4,
                culture: Culture {
                    binary_representation: 0b000_000_000_000_000,
                },

                effective_size: 5,
                number_offspring: 0,
                seasons_till_next_mutation: None,
                stored_resources: OneYearResources::from(1.),
                adaptation: ecology::Ecovector::from(p.minimum_adaptation),
            },
            Family {
                descendence: String::from("F"),
                location: start2,
                history: vec![],
                seasons_till_next_child: 4,
                culture: Culture {
                    binary_representation: 0b111_111_111_111_111,
                },

                effective_size: 5,
                number_offspring: 0,
                seasons_till_next_mutation: None,
                stored_resources: OneYearResources::from(1.),
                adaptation: ecology::Ecovector::from(p.minimum_adaptation),
            },
        ],
        t: 0,
        p: p,
    })
}

/**
# 6. Input Data

> Does the model use input from external sources such as data files or other
> models to represent processes that change over time?

 */
mod input {
    // None at the moment. It might come, though, when we use paleoclimate data.
}

/**
# 7. Submodels

> What, in detail, are the submodels that represent the processes listed in
> ‘Process overview and scheduling’? What are the model parameters, their
> dimensions, and reference values? How were submodels designed or chosen, and
> how were they parameterized and then tested?

*/
pub mod submodels {
    pub mod culture {
        use crate::Culture;
        use rand::prelude::*;

        pub fn distance(c1: Culture, c2: Culture) -> u32 {
            (c1.binary_representation ^ c2.binary_representation).count_ones()
        }

        pub fn mutate_culture(
            family_seasons_till_next_mutation: &mut Option<u32>,
            family_culture: &mut Culture,
            target_culture: Culture,
            culture_dimensionality: u8,
            culture_mutation_rate: f64,
        ) {
            *family_culture = target_culture;
            match family_seasons_till_next_mutation {
                None => {
                    *family_seasons_till_next_mutation =
                        Some(random::<f64>().log(1. - culture_mutation_rate) as u32);
                    mutate_culture(
                        family_seasons_till_next_mutation,
                        family_culture,
                        target_culture,
                        culture_dimensionality,
                        culture_mutation_rate,
                    );
                }
                Some(0) => {
                    let i: u8 = rand::thread_rng().gen_range(0, culture_dimensionality);
                    family_culture.binary_representation ^= 1 << i;
                    *family_seasons_till_next_mutation =
                        Some(random::<f64>().log(1. - culture_mutation_rate) as u32);
                }
                Some(k) => {
                    *family_seasons_till_next_mutation = Some(*k - 1);
                }
            }
        }
    }

    pub mod family_lifecycle {
        use crate::{adaptation, sensing};
        use crate::{Family, NodeId, OneYearResources, Parameters, Patch};
        use dashmap::DashMap;
        use rayon::prelude::*;

        pub fn procreate_and_migrate(
            families: &mut Vec<Family>,
            p: &Parameters,
            patches: &mut DashMap<NodeId, Patch>,
            cc: &DashMap<NodeId, usize>,
        ) -> Vec<Family> {
            families
                .par_iter_mut()
                .filter_map(|mut family| {
                    let nearby = sensing::nearby_locations(family.location, p);
                    // println!("{:?}", nearby);

                    let (destination, cost) = adaptation::decide_on_moving(
                        &family,
                        nearby.iter().map(|(i, d)| (*i, *d)),
                        patches,
                        false,
                        p,
                        cc,
                    );
                    // println!("Family {:} moved to {:}", family.descendence, destination);
                    
                    family.history.push(family.location);
                    family.location = destination;
                    family.stored_resources -= cost;

                    match maybe_procreate(&mut family, p.season_length_in_years) {
                        None => None,
                        // In terms of scheduling a new family can (and if possible
                        // should) move immediately when created. This behaviour is
                        // taken from del Castillo (2013).
                        Some(mut descendant) => {
                            let (destination, cost) = adaptation::decide_on_moving(
                                &descendant,
                                nearby.iter().filter_map(|(i, d)| {
                                    if *i == family.location {
                                        None
                                    } else {
                                        Some((*i, *d))
                                    }
                                }),
                                patches,
                                true,
                                p,
                                cc,
                            );
                            descendant.location = destination;
                            family.stored_resources -= cost;
                            Some(descendant)
                        }
                    }
                })
                .collect()
        }

        pub fn resources_at_season_end(
            resources: OneYearResources,
            size: usize,
            p: &Parameters,
        ) -> OneYearResources {
            resources - OneYearResources::from(size as f64) * p.season_length_in_years
        }

        pub fn maybe_procreate(family: &mut Family, season_length_in_years: f64) -> Option<Family> {
            if family.effective_size < 10 {
                None
            } else {
                family.number_offspring += 1;
                family.effective_size -= 2;
                Some(Family {
                    descendence: format!("{}:{:}", family.descendence, family.number_offspring),
                    location: family.location,
                    history: family.history.clone(),
                    seasons_till_next_child: (12. / season_length_in_years) as u32 + 1,
                    culture: family.culture,

                    effective_size: 2,
                    number_offspring: 0,
                    seasons_till_next_mutation: None,
                    stored_resources: OneYearResources::from(0.),
                    adaptation: family.adaptation,
                })
            }
        }

        pub fn maybe_grow(family: &mut Family, season_length_in_years: f64) {
            // print!("Growing {:} in {:}: size {:} ", family.descendence, family.seasons_till_next_child, family.effective_size);
            if family.seasons_till_next_child == 0 {
                family.effective_size += 1;
                family.seasons_till_next_child = (1.2 / season_length_in_years) as u32;
            } else {
                family.seasons_till_next_child -= 1;
            }
            // println!("becomes {:} ({:})", family.effective_size, family.seasons_till_next_child);
        }
        pub fn use_resources_and_maybe_shrink(
            size: &mut usize,
            resources: &mut OneYearResources,
            p: &Parameters,
        ) -> bool {
            let mut has_shrunk = false;
            while resources_at_season_end(*resources, *size, p) < OneYearResources::from(0.)
                && *size > 0
            {
                *size -= 1;
                has_shrunk = true;
            }
            *resources = resources_at_season_end(*resources, *size, p);
            has_shrunk
        }

        pub fn can_survive(family: &Family) -> bool {
            family.effective_size >= 2
        }
    }

    pub mod ecology {
        use crate::{OneYearResources, Parameters};

        /**

        This function models how much a family contributes to the group effort of
        harvesting a particular ecoregion. The family's contribution is given by its
        size, weighted with how well they know the ecoregion.

        As a side effect, because this function is the one that knows about the family's
        size, it also adds that size to a counter. This counter is used for the
        per-culture census. (The caller knows the culture of the family.)

        ```rust

        let mut family= model::Family::default();
        family.effective_size = 8;
        let mut sum_effort = 0;
        let (target, contribution) = model::raw_family_contribution(&mut family, 0, &mut sum_effort);
        assert_eq!(sum_effort, 8);
        *target += sum_effort as model::OneYearResources;
        assert_eq!(family.stored_resources, 8.0)
        ```
         */
        pub fn raw_family_contribution<'a>(
            family: &'a mut crate::Family,
            ecoregion: usize,
            sum_effort: &mut f64,
        ) -> (&'a mut OneYearResources, f64) {
            family.adaptation = family.adaptation * 0.95;
            family.adaptation.entries[ecoregion] += 0.05;
            let contribution = family.effective_size as f64 * family.adaptation[ecoregion];
            *sum_effort += contribution;
            (&mut family.stored_resources, contribution)
        }

        pub fn adjust_culture(
            group: &mut crate::collectives::Cooperative,
            p: &Parameters,
        ) -> crate::Culture {
            let target_culture = group.culture;
            for family in group.families.iter_mut() {
                crate::submodels::culture::mutate_culture(
                    &mut family.seasons_till_next_mutation,
                    &mut family.culture,
                    target_culture,
                    p.culture_dimensionality,
                    p.culture_mutation_rate,
                );
            }
            target_culture
        }

        pub fn group_contribution<'a>(
            group: &'a mut crate::collectives::Cooperative,
            ecoregion: usize,
            sum_effort: &mut f64,
        ) -> Vec<(&'a mut OneYearResources, f64)> {
            group
                .families
                .iter_mut()
                .map(|family| raw_family_contribution(family, ecoregion, sum_effort))
                .collect()
        }

        /**
         */
        pub fn exploit_patch<P>(
            mut groups: Vec<crate::collectives::Cooperative>,
            mut patch: P,
            p: &Parameters,
        ) where
            P: core::ops::DerefMut<Target = crate::Patch>,
        {
            for mut group in groups.iter_mut() {
                let total_storage: OneYearResources = group
                    .families
                    .iter_mut()
                    .map(|f| {
                        let s = f.stored_resources;
                        f.stored_resources = OneYearResources::from(0.0);
                        s
                    })
                    .sum();

                let mut sum_effort = 0.0;
                let mut contributions: Vec<_> = group_contribution(
                    &mut group,
                    crate::ecology::ATTESTED_ECOREGIONS,
                    &mut sum_effort,
                );

                for (storage, contribution) in contributions.iter_mut() {
                    **storage = total_storage * *contribution / sum_effort;
                }
            }

            for (ecoregion, (res, max_res)) in patch.resources.iter_mut() {
                let mut total_squared_groups_effort = 1.0;
                let mut work = groups
                    .iter_mut()
                    .map(|group| {
                        let mut sum_effort = 0.0;
                        let contributions: Vec<_> =
                            group_contribution(group, *ecoregion, &mut sum_effort);
                        total_squared_groups_effort += sum_effort.powi(2);
                        (sum_effort, contributions)
                    })
                    .collect::<Vec<(f64, Vec<(&mut OneYearResources, f64)>)>>();
                // println!("{:?} from {:?}", total_squared_groups_effort, work);
                work.drain(..)
                    .for_each(|(group_effort, mut contributions)| {
                        assert!(group_effort > 0.);
                        // println!("group effort: {:?}", group_effort);
                        let group_harvest =
                            *res * group_effort.powi(2) / total_squared_groups_effort;
                        assert!(group_harvest > OneYearResources::from(0.));
                        // println!("group harvest: {:?}", group_harvest);
                        let coefficient = std::cmp::min(
                            group_harvest / group_effort,
                            OneYearResources::from(2.5),
                        );
                        assert!(coefficient > OneYearResources::from(0.));
                        // println!("coefficient: {:?}", coefficient);
                        for (storage, contribution) in contributions.iter_mut() {
                            **storage += coefficient * p.season_length_in_years * *contribution;
                        }
                        let effort = coefficient * group_effort;
                        // println!("Harvest: {:?} (@ {:?}) from {:?}", effort, coefficient, res);
                        *res -= effort;
                        assert!(*res >= OneYearResources::from(0.));
                    });
                recover(res, *max_res, p.resource_recovery_per_season);
            }
        }

        /**
        Recover a patch, mutating its resources

        The recovery assumes exponential (or geometric, it's a discrete-time
                model) growth of resources in a patch up to a given maximum.

        The patch's `max_resouces` represents the maximum available for human
                extraction. This is a proportion `maximum_resources_one_adult_can_harvest` of the total
                primary production in that spot.

        For example, consider a patch with max_resources 2 kcal and current
                resources 1 kcal, where a quarter of the resources are accessible. This
                means that another 6 kcal is present, but inaccessible, for a total of
                7/8. The recovery is proportional to the current resources (i.e. 7 kcal)
                the ‘space left to grow’ (i.e. 1/8). With a resource recovery of 0.5,
                this means that this patch would recover by 7/16, so to 1.4375 available
                resources.

        ```rust
        # use model::submodels::ecology::recover;
        let mut resources = 1.0;
        recover(&mut resources, 2.0, 0.5, 0.25);
        assert_eq!(resources, 1.4375);
        ```

         */
        pub fn recover(
            patch_resources: &mut OneYearResources,
            patch_max_resources: OneYearResources,
            resource_recovery_per_season: f64,
        ) {
            if *patch_resources < patch_max_resources {
                *patch_resources = std::cmp::min(
                    patch_max_resources,
                    patch_max_resources * resource_recovery_per_season + *patch_resources,
                )
            }
        }
    }

    pub mod parameters {
        use serde_derive::{Deserialize, Serialize};
        use crate::movementgraph::MovementGraph;
        use crate::OneYearResources;

        #[derive(Debug, Serialize, Deserialize, Clone)]
        pub struct Parameters {
            pub resource_recovery_per_season: f64,
            pub culture_mutation_rate: f64,
            pub culture_dimensionality: u8,
            pub cooperation_threshold: u32,
            pub maximum_resources_one_adult_can_harvest: OneYearResources,
            pub evidence_needed: f64,
            pub payoff_std: f64,
            pub minimum_adaptation: f64,

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
                    maximum_resources_one_adult_can_harvest: OneYearResources::from(0.25),
                    evidence_needed: 0.1,
                    payoff_std: 0.1,
                    minimum_adaptation: 0.5,

                    season_length_in_years: 1. / 6.,
                    dispersal_graph: MovementGraph::default(),
                }
            }
        }
    }
}

pub fn store_state(state: State, statefile: String) -> Result<(), String> {
    let file = match File::create(statefile) {
        Ok(f) => f,
        Err(_) => return Err("Could not create state file".to_string()),
    };
    match serde_json::to_writer_pretty(file, &state) {
        Ok(_) => Ok(()),
        Err(_) => Err("Failed to store state.".to_string()),
    }
}

pub fn run(mut s: State, max_t: Seasons, o: &observation::ObservationSettings) {
    loop {
        step(&mut s.families, &mut s.patches, &s.p, s.t, o);

        if s.families.is_empty() {
            println!("Died out");
            break;
        }
        if (s.t == max_t) || (o.store_every > 0) && (s.t % o.store_every == 0) {
            println!("{:} % {:}", s.t, o.store_every);
            let to_be_stored = s.clone();
            let file = o.statefile.to_string();
            std::thread::spawn(|| {
                store_state(to_be_stored, file).unwrap_or_else(|r| println!("{:}", r))
            });
        }
        if s.t >= max_t {
            println!("Ended");
            break;
        }
        s.t += 1;
    }
}
