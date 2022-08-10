/*!

Model Description
=================

This model description follows the ODD (Overview, Design concept, Details)
protocol (Grimm et al., 2006; Grimm et al., 2010). Where appropriate, we quote
the guiding questions of the ODD protocol at the beginning of a section. The
model description follows the idea of literate programming (Knuth 1992) to the
extent useful in Rust source code – the actual model is generated from this file
that documents the model, but the source code is largely commented using
natural-language descriptions, not generated from them (as would be the case in
a literal program).

## 1. Purpose

The dispersal model generates a deep phylogeny of hunter-gatherer cultures based
on culture-mediated cooperation and resource-driven migration. It is a
demographic migration model in which the areal distribution of languages is an
emergent property, not an imposed structure. As such, there is a feedback loop
between resource availability, linguistic proximity, and cooperation, which
drives the population dynamics of hunter-gatherers in the model.

The model is designed with extension to more concrete research questions in
mind. In the current, first stage, the purpose of the model is to investigate
how languages disperse and split, driven only by the necessary interactions
between humans. The particular research question is whether distinct cultural
areas, with reasonably sharp borders, can emerge from only a simple feedback
loop between conditional cooperation, linguistic assimilation, and genetic
drift.

If that is the case, the summary statistics of a phylogeny (in particular
diversification rates) can later be compared to values known from language
evolution. The model is structured to be easily applied to study the history of
the settlement of the Americas at a later time. It would require paleoclimate
data and interpretations of past ecoregions (and potentially an extension
describing the changes in human behaviour patterns) to produce results that can
be compared to that history.

To construct a demographic migration model with culture, for the purpose of the
project we set out here, we need the following ingedients.

1. A representation of culture, which is at the least able to undergo neutral evolution
   (showing heritability and random variation, not necessarily fitness) and which in
   addition allows horizontal transfer of cultural traits other than from a parent to a
   child population.
2. Agents that carry cultural traits and are located in geographical space, which has
   differing ecological features
3. A system that drives the demographics of the agents in time and space
4. A way for culture and population dynamics to interact in a way that can create create
   distinct cultural areas instead of a vast cultural cline or dialect continuum.
   In the following subsections, we will consider each of these elements separately.

*/

// TODO: Dear Rustacean, I know that my use of documentation comments is
// hazardous, because they land in the generated documentation in a different
// order, or sometimes attached to things that are only accidentally coming afterwards.

// Load useful modules

#![allow(clippy::redundant_field_names, clippy::implicit_hasher)]

use dashmap::DashMap;
use serde_derive::{Deserialize, Serialize};

use rustc_hash::{FxHashMap, FxHasher};
use std::fs::File;

use rayon::prelude::*;
use rayon::iter::{IntoParallelIterator, IntoParallelRefMutIterator};

pub mod argparse;
mod debug;
pub mod ecology;

use submodels::parameters::Parameters;

/**

## 2. Entities, state variables, and scales

The model consists of agents interacting on a network of habitable patches in
discrete time. One time step is supposed to model a season during which a group
of hunter-gatherers remains within the same area. Some of the underlying data is
using seconds as base unit and needs to be translated.

 */
pub type Seasons = u32;
const SECONDS_PER_YEAR: f64 = 365.24219 * 24. * 60. * 60.;

/**
Whereever possible, resources are measured in terms of the resources one adult
consumes in a year, abbreviated OYR. For comparison, in terms of food energy,
one OYR corresponds to about 650'000 kcal/3'000'000 kJ, but it also includes
other limited resources necessary for human survival.

 */
use ecology::OneYearResources;

/**
### 2.1 Geography: Nodes and Patches

The geography of the simulation is described by a directed weighted graph, where
each edge has a weight representing the travel time (by foot, or by simple boat
along rivers and coasts) in seconds. This movement graph is constructed before
the initialization from pre-processed real geographical data. The construction
of our movement graph is described in [@kaiping2021network]. Each node has a
unique numerical identifier.

 */
pub mod movementgraph;
pub type NodeId = petgraph::graph::NodeIndex<usize>;
/**

It also contains the data on the patch of land it corresponds to. A patch
contains one or more terrestrial ecoregions. For each ecoregion, the node tracks
maximum and current availability of resources, measured in OYR available over
the course of one season. We estimate these numbers from population densities,
adapting the estimation from [@tallavaara2018productivity] to derive population counts, and
therefore the resource amounts necessary to sustain them, in the separate
ecoregions areas that are most easily reached from the point represented by the
node.

Our areas have an average area of about 400 km² each, but because we inlude all
locations best accessibility from each hexagon center (in a generalization of
Voronoy polygons, see [@kaiping2021network]) instead of the hexagons as defined
in the underlying regular grid, the actual area represented by a patch can vary:
On the boundary between a steep, inaccessible canyon and a flat grassland the
patches inside the canyon may be much smaller than the patches on the grassland.

In the current state of the model, the patch at a node is constant throughout
the simulation, but future work might add seasonality or random variation. By
nature, the simulation invites the extension to include paleoclimate data, such
as recently published in [@beyer2020highresolution], but that only becomes
relevant once it shows fundamentally reasonable dynamics.

 */
#[derive(Serialize, Deserialize, Clone)]
pub struct Patch {
    /// For every relevant terrestrial ecoregion area accessible from the node
    /// location, the tuple of resources available for the next season, and the
    /// maximum possible resources available in a season, both in OYR.
    resources: FxHashMap<usize, (OneYearResources, OneYearResources)>,
}

/**
### 2.2 Families

The main decision-making agents of the simulation are families
[@delcastillo2013modeling;@barcelo2014social]. Families can migrate between
nodes, grow, or shrink. New families can split off from existing families and
become independent agents, and families that shrink below 2 individuals die out
and are remove from the simulation. Each family has the following properties.

> MAYBE: Lake (2000)

 */

#[derive(Serialize, Deserialize, Clone)]
pub struct Family {
/**
- A history of decendence.

     */
    pub descendence: String,
    /**
- The effective size of the family in number of adults. One adult is assumed to
    consume the same amount of food/energy and to contribute the same labor to
    foraging as any other adult. For simplicity, children are not modelled, even
    though their needs are known to have an influence on human cultural patterns
    [@surovell2000early;@volk2013infant].

     */
    pub effective_size: usize,
    /**
- A culture shared among all members of the family, to be detailed below

     */
    pub culture: Culture,
    /**
- Its current location, the numerical ID of a node in the movement graph.

     */
    pub location: NodeId,
    /**
- The amount of resources, in OYR per capita, which the family harvested in the previous season.

     */
    pub last_harvest: OneYearResources,
    /**
- The amount of stored resources, in OYR, the family has access to without going foraging.

     */
    pub stored_resources: OneYearResources,
    /**
- Its adaptation to local ecoregions is represented as a vector with values
      between 0.0 (completely unknown) and 1.0 (perfectly familiar) for any
      ecoregion the family might encounter.

     */
    pub adaptation: ecology::Ecovector,
    /**
- The number of seasons to wait until the next (unsimulated) child becomes an adult.

     */
    pub seasons_till_next_adult: Seasons,
    /**

- For bookkeeping purposes (eg. generating descendant's ‘descendence’ values),
      keep track of the number of families that have split off from this family so
      far.

     */
    pub number_offspring: u16,
    /**

- The previous locations of the agent. The most recent locations serve as the memory of the agent. It is also useful for some bits of analysis

     */
    pub history: Vec<(NodeId, OneYearResources)>,
    /**
- A bookkeeping quantity. Instead of mutation happening in each time step with
      a tiny probability, the same distribution of times between mutations is
      generated by drawing a number of season from a geometric distribution and
      counting it down by one each time step, which is useful because random number
      generation is computationally expensive.

     */
    pub seasons_till_next_mutation: Option<Seasons>,
}
/**

The descendence of a family also serves as its unique identifier.

 */
impl PartialEq for Family {
    fn eq(&self, other: &Self) -> bool {
        self.descendence == other.descendence
    }
}

/**
### 2.3 Bands

Families that have compatible cultures and are in the same location can band
together to share resources between different families, and to compete against
other bands. The cooperative bands formed by cooperating families are
higher-level agents created ad-hoc in each time step. They do not persist or
have effect beyond a single time step. Resource exploitation happens at the
level of the band and is distributed to the individual families after the fact.

 */
pub struct Band<'a> {
    pub families: Vec<&'a mut Family>,
    pub culture: Culture,
}
/**
### 2.4 Cultures

Every family has a culture. These are very abstract and vastly simplified, due
to the lack of quantitative data on cultural evolution in a framework comparable
to the one used for this study. Based on the need to have a data structure that
supports random drift, a low but positive chance of back-mutation, and a large
number of equally (dis-)similar cultures, we describe culture using a binary
vector. Similar models of culture have been used in [@debie2007agentbased]. A more
detailed discussion of the choices and interactions can be found in [Submodel:
Culture].

The cultures can be observed, so the implementation defines binary and
hexadecimal representations for cultures.

 */
#[derive(Ord, PartialOrd, Eq, PartialEq, Hash, Clone, Serialize, Deserialize)]
pub struct Culture {
    // TODO: Once integer generics are a thing in Rust, change from Vec to an
    // array of a defined constant size.
    in_memory: Vec<usize>,
}

impl std::fmt::Binary for Culture {
    fn fmt<'a>(&self, f: &mut std::fmt::Formatter<'a>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            self.in_memory
                .iter()
                .map(|i| format!("{:b}", i))
                .collect::<Vec<String>>()
                .join("")
        )
    }
}

impl std::fmt::LowerHex for Culture {
    fn fmt<'a>(&self, f: &mut std::fmt::Formatter<'a>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            self.in_memory
                .iter()
                .map(|i| format!("{:x}", i))
                .collect::<Vec<String>>()
                .join("")
        )
    }
}

/**
### 2.5 State

The status of all families, including their locations, is the core of the model
state. The state also tracks time, measured in time steps corresponding to one
season each, since the start of the simulation, and stores the mapping between
nodes and patches as well as a copy of the model parameters.

 */
#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct State {
/**

- The agents (families) currently active in the simulation. This vector
    changes over time as new families are added and old families die out.

     */
    families: Vec<Family>,
    /**

- The current time step, in half years since start of the simulation.

     */
    pub t: Seasons,
    /**

- The patches of the model, indexed by node ID in a graph. This set is
    fixed, although the properties of the patches may change with time.

     */
    patches: DashMap<NodeId, Patch>,
    /**

- The parameters of the model.

     */
    p: Parameters,
}
/**

## 3. Process overview and scheduling

> Who (i.e., what entity) does what, and in what order? When are state
> variables updated? How is time modeled, as discrete steps or as a continuum
> over which both continuous processes and discrete events can occur? Except
> for very simple schedules, one should use pseudo-code to describe the
> schedule in every detail, so that the model can be re-implemented from this
> code. Ideally, the pseudo-code corresponds fully to the actual code used in
> the program implementing the ABM. [@grimm2010odd]

The model progresses in discrete time steps, each corresponding to a fixed fraction of a year.
The structure of a single time step consist of two parts, as follows.

 */
fn step(
    families: &mut Vec<Family>,
    patches: &mut DashMap<NodeId, Patch>,
    storage: &DashMap<
        NodeId,
        FxHashMap<Culture, (usize, OneYearResources)>,
        std::hash::BuildHasherDefault<FxHasher>,
    >,
    p: &Parameters,
    t: Seasons,
    nearby_cache: &mut DashMap<NodeId, Vec<NodeId>, std::hash::BuildHasherDefault<FxHasher>>,
    o: &observation::Settings,
) -> DashMap<
    NodeId,
    FxHashMap<Culture, (usize, OneYearResources)>,
    std::hash::BuildHasherDefault<FxHasher>,
> {
    update_each_family_step(families, patches, storage, nearby_cache, p);
    distribute_each_patch_resources_step(families, patches, p, o, t)
}
/**
The first part focuses on the individual families, which shrink, grow, die,
split, and move. It constructs the mapping of families at the end of the season,
grouped by their location at the end of the movement. Because the movement of a
family depends on the distribution of the families at he start of the season and
not at the time of movement, it can happen entirely in parallel.

 */
fn update_each_family_step(
    families: &mut Vec<Family>,
    patches: &mut DashMap<NodeId, Patch>,
    storage: &DashMap<
        NodeId,
        FxHashMap<Culture, (usize, OneYearResources)>,
        std::hash::BuildHasherDefault<FxHasher>,
    >,
    nearby_cache: &mut DashMap<NodeId, Vec<NodeId>, std::hash::BuildHasherDefault<FxHasher>>,
    p: &Parameters,
) {
    let cc: DashMap<NodeId, usize, std::hash::BuildHasherDefault<FxHasher>> = DashMap::default();
    for f in families.iter() {
        *cc.entry(f.location).or_insert(0) += f.effective_size;
    }

    families.par_iter_mut().for_each(|mut family| {
        if submodels::family_lifecycle::use_resources_and_maybe_shrink(
            &mut family.effective_size,
            &mut family.stored_resources,
            p,
        ) {
            // TODO: Think some more about what happens to the children logic when you starve
            // For now, for simplicity's sake, all children starve.
            family.seasons_till_next_adult = (15. / p.season_length_in_years) as u32 + 1;
        }

        submodels::family_lifecycle::maybe_grow(family, p.season_length_in_years);
    });

    families.retain(submodels::family_lifecycle::can_survive);

    let mut children = submodels::family_lifecycle::procreate_and_migrate(
        families,
        storage,
        p,
        patches,
        &cc,
        nearby_cache,
    );
    families.append(&mut children);
}

// IDEAS:
// every agent stores past interactions, plays tit-for-tat (C→D→L→C)
// Start at “cooperate with self”
// (choose to play cooperator/defector/loner in public goods game to extract resources)
// Store a base strategy, starting L
// Base strategy changes HOW? Random mutation? Or maybe to some function of history?
/**
The second part focusses on the patches. The resources of a patch are updated
according to the families exploiting them over the season. This is described in
detail in [Submodel Exploitation]. Everything here happens locally to a patch,
with no external interaction, so this can be done in parallel. After
exploitation, patches recover for the next season according to [Submodel Patch].
This concludes a time step.

 */
fn distribute_each_patch_resources_step(
    families: &mut Vec<Family>,
    patches: &mut DashMap<NodeId, Patch>,
    p: &Parameters,
    o: &observation::Settings,
    t: Seasons,
) -> DashMap<
    NodeId,
    FxHashMap<Culture, (usize, OneYearResources)>,
    std::hash::BuildHasherDefault<FxHasher>,
> {
    let families_by_location: DashMap<
        NodeId,
        Vec<&mut Family>,
        std::hash::BuildHasherDefault<FxHasher>,
    > = DashMap::default();
    families.par_iter_mut().for_each(|f| {
        families_by_location
            .entry(f.location)
            .or_insert_with(Vec::default)
            .push(f)
    });

    let stored_resources = DashMap::default();
    let cultures_by_location: FxHashMap<NodeId, FxHashMap<Culture, usize>> = families_by_location
        .into_par_iter()
        .map(|(patch_id, families)| {
            let mut cc = FxHashMap::default();

            let mut groups = crate::collectives::cooperate_or_fight(families, p);

            let mut cs = stored_resources
                .entry(patch_id)
                .or_insert_with(FxHashMap::default);
            for group in groups.iter_mut() {
                submodels::ecology::adjust_culture(group, p);
                let group_size = group.families.iter().map(|f| f.effective_size).sum();
                cc.insert(group.culture.clone(), group_size);
                cs.insert(
                    group.culture.clone(),
                    (group_size, OneYearResources::default()),
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
    stored_resources
}

/**
## 4. Design concepts

Under the ODD protocol, the design principles largely fall into questions. Where
my answer indicates an invariant property of the simulation, I provide a test
function that checks that invariance if possible.

 */
mod concepts {

/**
### 4.1 Basic priciples

> Which general concepts, theories, hypotheses, or modeling approaches are
> underlying the model’s design?

> How were they taken into account?

> Are they used at the level of submodels, or is their scope the system level?

> Will the model provide insights about the basic principles themselves, i.e.,
> their scope, their usefulness in real-world scenarios, validation, or
> modification?

> Does the model use new, or previously developed, theory for agent traits from
> which system dynamics emerge?

*/
    mod basic_principles {
/**

To achieve realistic outcomes comparable, at least in part, to real human
history, the model presented here is not the most simple abstract model
including the components listed below, but attempts to represent these submodels
realistically in as far as existing data (and necessary simplifications to
reduce model complexity) permit.
    
Fundamentally, the model presented here is a demographic dispersal model. The
minimal conditions on the model given by this foundation are the following two:
1. Agents must be able to **move through the continental geography**
(potentially at a cost) 2. The model includes **local population capacities**,
which may be realized as soft bounds which permit the actual population size to
overshoot the carrying capacity temporarily.

Population dynamics and migratory patterns are to a large extent driven by the
local geographical conditions. Difficult terrain impedes movement according to
the imposed graph structure, which has edges with a higher cost where
cross-terrain movement is more difficult, and marginal environments are less
attractive to migrating agents and depleted quicker.

 */
    /* TODO: Write a test where a single patch has a given carrying capacity,
     * and a cooperating population grows to that capacity. */
/**

According to its purpose of providing a model for language dispersal and split,
our model draws on existing publications looking at cultures changing in time
and space, focussing on taking a bottom-up approach on languages splitting. A
major useful reference is the PSMED model [@delcastillo2013modeling;
@barcelo2013psmed; @barcelo2015simulating], which includes two drivers of
cultural change in the agents: Drift and assimilation. In isolation, agents'
**cultures undergo drift**. Where social interactions happen, drift is
counterbalanced by cultural assimilation for agents that cooperate with each
other. In the case of ths PSMED model, agents interact with other agents in a
limited range and based on need. We posit a similar basic relationship between
**cooperation driving cultural assimilation, and cultural similarity
conditioning cooperation**. Due to the larger scope of the model presented
here, however, we avoid tracking networks of pairwise interactions between
agents. We instead assume that cooperation happens within patches, where
groups of families in the patch form based on agents' cooperation strategies.
These groups participate in a **public goods game** [@kagel2020handbook],
allowing for cooperation and the sharing of resources. Instead of the network
between agents, our families track cultures they are willing to cooperate
with. Our maximum range of interactions is about 20 km due to the size of our
patches.

A major necessary feature of the dispersal component for the culture component
is to allow **group fission and fusion**, to allow culture splits to emerge from
low-level agent interactions, but to also drive incentives for agents to
congregate, cooperate, and assimilate their cultures. As such, the model draws
heavily on the agent-based model by [@crema2014simulation; @crema2015modelling].

One of the deviations from Crema's model, and also from PSMED, is the geography
underlying the simulation. Crema and Barceló use a regular quadratic grid with
arbitrary fixed or randomly generated resources. A long-term goal for our model
it applicability to real-world language dispersal processes. As such, we derive
the resources from real-world data and use a **network derived from a hexagonal
grid** to inform the topology. Such a grid has some advantages over a quadratic
grid on the plane, in that all hexagons within a fixed number of steps from a
focus are a much better approximation of a circle than all squares within a
fixed number of steps from a central square. Hexagons also have nicer properties
for delineating boundaries, which need some approximation for square grids (cf.
[@kuijper2004detecting], which is not directly on finding boundaries, but
solving that problem nonetheless). On a sphere such as Earth, using a hexagonal
grid has the additional advantage of reducing distortions. Using hexagonal grids
in continental or global dispersal simulations is well-established
[@gavin2017processbased; @callegari2013agentbased].

Our model is thus an agent-based model of the emergence of cultural areas from
fission-fusion dynamics, where local cooperation is the main driver to prevent
cultures from unconstrained evolutionary drift.

It has been hypothesized that there are different demographic dispersal patterns
for different periods in the settlement of the Americas. Presumably, the first
humans dispersing in the Americas were spreading very fast, and various
mechanisms have been proposed to explain this (eg. leapfrog spread to
high-quality patches [@prates2020rapid;citations-therein], technology
specializing on hunting migratory big game [jennings2008san], or a “kelp
highway” facilitating southward movement along the west coast of the Americas
[@erlandson2007kelp;@erlandson2015kelp]). For the early pre-historic
hunter-gatherer groups, when their immediate neighborhood is persistently
settled, basic interaction strategies appear to be important []. Where groups
are somewhat stable, where finding suitable sexual partners becomes a major
constraint, with different strategies necessary in marginal regions (where the
limiting factor is the opportunities to meet potential partners []) and in more
densely populated regions (where linguistic exogamy develops, partly to avoid
inbreeding []). Lastly, the “farming/dispersal hypothesis”
[robbeets2017language] suggests that the development of agriculture may have led
to the spread of some language families, replacing others, in conjunction with
an increased population capacity.

Our simulation focusses on the second stage, where the continents are **largely
settled, but the local structure still fluctuates**. We take into account some
aspects of the leapfrog theory of settlement, by assuming that **resources
between patches may vary greatly**, and where multiple implementation options
are available, we choose the option that is the better approximation of sexual
partner choice. We **disregard farming** and the associated changes in
population capacity. We also ignore the different ecological and climatic state
of the uninhabited Americas which would be necessary to investigate different
settlement time depsths and to compare the different early dispersal models.

 */
    }
/**

### 4.2 Emergence

> What key results or outputs of the model are modeled as emerging from the
> adaptive traits, or behaviors, of individuals? In other words, what model
> results are expected to vary in complex and perhaps unpredictable ways when
> particular characteristics of individuals or their environment change?

*/
    mod emergence {
        use crate::*;
/**

The main emergent property we aim for will be culturally somewhat uniform
territories. We expect that the interplay of migration, cooperation and cultural
similarity leads to regions sharing similar cultures, with noticeable
boundaries. This means that the plot of cultural distances vs. geographical
distances should therefore show small cultural distances for small geographical
distances. There should be a critical geographical distance of cohesion where
the distribution of cultural distances becomes bimodal, with one mode being
lower than the cooperation threshold and one mode above the cooperation
threshold. For large geographical distances, the cultural distances should be
high, but with a big variance. The critical geographical distance is likely to
depend on the region, being larger in more marginal environments where migration
is more frequent.

Underlying the model is a fission-fusiom model of group dynamics
[@crema2014simulation], so we expect a similar analysis to apply to our model.
Crema observes actual fission-fusion cycles only for a small part of the
parameter space, so it is not clear whether we should expect them for our model.

Cooperation in the model is part of the strategy choices available to the
agents. It is an important question in theoretical biology and related fields
under what circumstances cooperation can prevail over defectors and other
free-riders, and it may even affect the dispersal of languages (though more
likely on the level of collectives, where the interactions between groups of
different cultures range from assimilation and language shift all the way to war
and lethal violence). One regime sustaining cooperation seems to be the
inclusion of “Loner” strategies in Public Goods Games (PGG) [@?], and because
this strategy corresponds to the avoidance of other agents and thus potentially
to strict boundaries, we include a PGG in the model structure [Submodel ??] and
expect varying levels of cooperation to emerge in the model.

*/

    /** Compute the Hamming distance between two culture vectors */
    #[cfg_attr(target_arch = "x86_64", target_feature(enable = "popcnt"))]
    pub unsafe fn distance(c1: &Culture, c2: &Culture) -> u32 {
        c1.in_memory
            .iter()
            .zip(c2.in_memory.iter())
            .map(|(bits1, bits2)| (bits1 ^ bits2).count_ones())
            .sum()
    }

    pub fn will_cooperate(c1: &Culture, c2: &Culture, cooperation_threshold: u32) -> bool {
        unsafe { distance(c1, c2) < cooperation_threshold }
    }
}
}

/**
> Are there other results that are more tightly imposed by model rules and hence
> less dependent on what individuals do, and hence ‘built in’ rather than
> emergent results?

Adaptation to local environment is implemented to grow with exposure time, without underlying mechanisms being reflected, as described below.

### 4.3 Adaptation

 */
mod adaptation {
    use crate::*;
/**

> What adaptive traits do the individuals have? What rules do they have for
> making decisions or changing behavior in response to changes in themselves or
> their environment? Do these traits explicitly seek to increase some measure of
> individual success regarding its objectives?

Agents have control over two main traits. The first of these is their location.
Agents optimize (within the limits of their knowledge) their location to
increase the amount of resources available for them to gather. The resources
gathered indirectly contribute to future success: A certain minimum of gathered
resources is necessary to procreate, but extremely large amounts gathered soon
do not increase the immediate rate of procreation and can in extreme cases
(overexploitation of local resources) even be detrimental in the long run and
through sharing instead benefit other families of the same culture.

 */
    use std::sync::Arc;

    pub fn decide_on_moving<'a>(
        family: &Family,
        storage: &DashMap<
                NodeId,
            FxHashMap<Culture, (usize, OneYearResources)>,
            std::hash::BuildHasherDefault<FxHasher>,
            >,
        known_destinations: &'a [NodeId],
        patches: &'a DashMap<NodeId, Patch>,
        quality_cache: &Arc<DashMap<NodeId, OneYearResources>>,
        p: &Parameters,
        population: &DashMap<NodeId, usize, std::hash::BuildHasherDefault<FxHasher>>,
    ) -> (NodeId, OneYearResources) {
        let travel_costs = sensing::many_nearby_locations(family.location, p);

        // 0.5 = memory_decay_per_season ^ memory_half_life
        // Assuming a half life of 2 years:
        let memory_decay_per_season = (0.5_f64).powf(p.season_length_in_years / 2.0);
        let mut memory = 1.0;

        let destination_expectation =
            known_destinations
            .iter()
            .map(|destination_id| (destination_id, None))
            .chain(family.history.iter().rev().filter_map(
                |(destination_id, expected_harvest)| {
                    if memory < f64::EPSILON {
                        return None;
                    }
                    if fastrand::f64() < memory {
                        memory *= memory_decay_per_season;
                        Some((destination_id, Some(expected_harvest)))
                    } else {
                        memory *= memory_decay_per_season;
                        None
                    }
                },
            ))
            .filter_map(|(destination_id, expected_harvest)| {
                if let Some(cost) = travel_costs.get(destination_id) {
                    let harvest = if let Some(x) = expected_harvest {
                        *x
                    } else {
                        objectives::destination_quality(
                            destination_id,
                            family,
                            storage,
                            *population.get(destination_id).as_deref().unwrap_or(&0),
                            patches.get(destination_id).as_deref().unwrap(),
                            p,
                            quality_cache,
                        )
                    };
                    Some((destination_id, harvest, cost))
                } else {
                    None
                }
            });
        best_location(destination_expectation, family.location)
    }
/**
The other adaptive trait of a family is their strategy for interacting with
other families. Every family has a basic strategy to be used with unknown
cultures, as well as a fixed memory storing how they intend to interact with
families they have encountered before. These strategies are modified as follows.

 - When a family ...

 */

    pub fn update_interaction_strategy() {}

/**

In order to account for limited knowledge in an economic fashion without
tracking the kowledge of every single family, and to avoid the case where one
location is the best choice of a large amount of agents by a very small margin,
agents choose the ‘best’ location (according to expected resources with a
discount for enemies, as explained below) subject to a random re-scaling of the
expected quality by a factor uniformly distributed between 0 and 1.

*/
    pub fn best_location<'a, KD>(kd: KD, null: NodeId) -> (NodeId, OneYearResources)
    where
        KD: IntoIterator<Item = (&'a NodeId, OneYearResources, &'a OneYearResources)>,
    {
        let (mut target, mut max_gain, mut t_cost) = (
            &null,
            OneYearResources::from(0.0),
            OneYearResources::from(0.0),
        );
        let mut m = max_gain - t_cost;

        for (location, gain, l_cost) in kd {
            // println!("Location {:?} would provide {:?} but cost {:?} to move to.", location, gain, l_cost);
            let g = gain * stochasticity::resource_uncertainty() - *l_cost;
            if g > m {
                target = location;
                max_gain = max_gain.max(gain);
                t_cost = *l_cost;
                m = max_gain - t_cost;
            }
        }
        (*target, t_cost)
    }
/**

> Or do they instead simply cause individuals to reproduce observed behaviors
> that are implicitly assumed to indirectly convey success or fitness? [@grimm2010odd]

It has been proposed [@XXX, maybe-somewhere-bantu] that the need to adaptat to
new environments slows migrations on a continental scale. As a consequence, each
family has a vector listing their success at extracting resources from the
various ecoregions. This is not modeled as an evolutionary adaptive process;
instead, the environmental adaptation of a family goes up deterministically when
they interact with an ecoregion, and goes down to a (parameterized) minimum
otherwise. This submodel is described in [7.X: Environmental adaptation].

*/
}

/**
### 4.4 Objectives

> If adaptive traits explicitly act to increase some measure of the individual’s
> success at meeting some objective, what exactly is that objective and how is
> it measured?
>
> When individuals make decisions by ranking alternatives, what criteria do they
> use? [@grimm2010odd]

 */
mod objectives {
    use crate::*;
/**

Agents base their decision for one potential destination of migration over
another based on a modification of the resources expected per individual at the
destination.

The amount of resources expected per individual is the sum of

 - the resouces stored by the ‘friendly’ families (i.e. families with compatible
    language) currently in the region, divided by the number of ‘friends’ and

    */
    pub fn from_friends(
        i: &NodeId,
        family: &Family,
        storage: &DashMap<
            NodeId,
            FxHashMap<Culture, (usize, OneYearResources)>,
            std::hash::BuildHasherDefault<FxHasher>,
        >,
        p: &Parameters,
    ) -> (usize, OneYearResources) {
        let (friends, stored_by_friendlies) = match storage.get(i) {
            None => (0, OneYearResources::from(0.0)),
            Some(h) => {
                match h.iter().find(|c| {
                    emergence::will_cooperate(&family.culture, c.0, p.cooperation_threshold)
                }) {
                    None => (0, OneYearResources::from(0.0)),
                    Some((_, v)) => *v,
                }
            }
        };

        let stored_by_friendlies_except_me = if *i == family.location {
            stored_by_friendlies - family.stored_resources
        } else {
            stored_by_friendlies
        };
        (
            friends
                + (if family.location == *i {
                    0
                } else {
                    family.effective_size
                }),
            stored_by_friendlies_except_me,
        )
    }
    /**

 - for each ecoregion in the patch, the resources accessible to the family this
       and next season, divided by the total population of the patch at the
       beginning of the time step. Note that here the family extrapolates its own
       environmental adaptation to the other families in the patch.

    */
    pub fn expected_harvest(patch: &Patch, p: &Parameters) -> OneYearResources {
        let harvestable_resources = patch.resources.values();
        harvestable_resources
            .map(|(res, max_res)| (*res + *max_res * p.resource_recovery_per_season))
            .sum::<OneYearResources>()
    }

    pub fn has_enemies(
        i: &NodeId,
        culture: &Culture,
        storage: &DashMap<
            NodeId,
            FxHashMap<Culture, (usize, OneYearResources)>,
            std::hash::BuildHasherDefault<FxHasher>,
        >,
        cooperation_threshold: u32,
    ) -> bool {
        match storage.get(i) {
            None => false,
            Some(population) => {
                for other_culture in population.keys() {
                    if !emergence::will_cooperate(culture, other_culture, cooperation_threshold) {
                        return true;
                    }
                }
                false
            }
        }
    }

    use std::sync::Arc;
    /**

This expected quality is modified to account for the risk of encountering
enemies. If the cooperative families are not the majority in the destination,
there is a significant danger from warfare to the migrating family. They adjusts
the expected quality of the destination by a constant factor (a model parameter)
for every enemy that is in excess of the family's language being the dominant
one in the region.

*/
// A family would usually have a size around 5, I guess, and each family
// encountered has probability 1/2 to mess with this family. So each individual
// has a chance of 0.5^0.2 – and that's the default of the parameter.

    pub fn destination_quality<'a>(
        i: &'a NodeId,
        family: &'a Family,
        storage: &DashMap<
            NodeId,
            FxHashMap<Culture, (usize, OneYearResources)>,
            std::hash::BuildHasherDefault<FxHasher>,
        >,
        population: usize,
        patch: &Patch,
        p: &Parameters,
        quality_cache: &Arc<DashMap<NodeId, OneYearResources>>,
    ) -> OneYearResources {
        let pop_at_destination = population
            + (if family.location == *i {
                0
            } else {
                family.effective_size
            });

        if has_enemies(i, &family.culture, storage, p.cooperation_threshold) {
            OneYearResources::default()
        } else {
            *quality_cache
                .entry(*i)
                .or_insert(expected_harvest(patch, p))
                / pop_at_destination as f64
        }
    }
}
/**
### 4.5 Learning

> M[ay] individuals or agents (but also organizations and institutions) change
> their adaptive traits over time as a consequence of their experience? If so,
> how? [@grimm2010odd]

INTERACTION STRATEGY UPDATING
 */

mod learning {}

/**
### 4.6 Prediction

> Prediction is fundamental to successful decision-making; if an agent’s
> adaptive traits or learning procedures are based on estimating future
> consequences of decisions, how do agents predict the future conditions (either
> environmental or internal) they will experience? If appropriate, what internal
> models are agents assumed to use to estimate future conditions or consequences
> of their decisions? What tacit or hidden predictions are implied in these
> internal model assumptions? [@grimm2010odd]

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
### 4.7 Sensing

> What internal and environmental state variables are individuals assumed to
> sense and consider in their decisions? What state variables of which other
> individuals and entities can an individual perceive; for example, signals that
> another individual may intentionally or unintentionally send? Sensing is often
> assumed to be local, but can happen through networks or can even be assumed to
> be global (e.g., a forager on one site sensing the resource levels of all
> other sites it could move to). If agents sense each other through social
> networks, is the structure of the network imposed or emergent? Are the
> mechanisms by which agents obtain information modeled explicitly, or are
> individuals simply assumed to know these variables? [@grimm2010odd]

 */
mod sensing {
    use crate::{movementgraph, NodeId, OneYearResources, Parameters, SECONDS_PER_YEAR};
/**

Individuals know about nearby locations. The exploration of those locations is
not explicitly modelled.

*/
// I should say a bit more about sensing, I think.
    use std::collections::HashMap;
    pub fn many_nearby_locations(
        location: NodeId,
        p: &Parameters,
    ) -> HashMap<NodeId, OneYearResources> {
        // Triple: opportunity cost from not foraging, and cost from doing likely
        // heavy labor instead.
        let cost_per_walking_second: OneYearResources =
            OneYearResources::from(1.) / SECONDS_PER_YEAR * 3.;

        movementgraph::bounded_dijkstra(&p.dispersal_graph, location, MAX_SEASON_RANGE, |e| {
            *petgraph::graph::EdgeReference::weight(&e)
        })
        .drain()
        .filter_map(|(n, d)| {
            if d.is_normal() {
                Some((n, cost_per_walking_second * d))
            } else {
                None
            }
        })
        .collect()
    }

    use petgraph::visit::EdgeRef;
    use std::collections::HashSet;
    pub fn nearby_locations(location: NodeId, p: &Parameters) -> Vec<NodeId> {
        p.dispersal_graph
            .edges(location)
            .flat_map(|e| {
                let neighbor = if e.source() == location {
                    e.target()
                } else {
                    e.source()
                };
                p.dispersal_graph.edges(neighbor).map(move |ne| {
                    if ne.source() == neighbor {
                        ne.target()
                    } else {
                        ne.source()
                    }
                })
            })
            .collect::<HashSet<NodeId>>()
            .drain()
            .collect()
    }
/**

#### Technical details on the maximum sensing distance

There is a maximum range on sensing: [@kelly2013, Table 4.1] (reproduced in the
supplementary material) lists mobility data for 70 hunter-gatherer groups,
including figueres for the yearly territory of a group where that number could
be inferred. Three quarters of the listed territories are about 2000 km² or
less. In our model, this corresponds to roughly 5 patches (each with an area of
ca. 400 km²) being visited per year. To give some leeway, we permit our agents
up to 6 moves per year, so a season corresponds to two months. The third
quartile, across all groups, of the yearly migratory distance is around 320 km,
so in this simulation the maximum distance moved per season is about 53.3 km.
The simulation is centered around distances measured in time units (seconds,
seasons, years). In all biomes except for mangroves, the terrain speed factor is
greater than one (see `supplement/distances/ecoreginos.py`). If humans generally
prefer 2% downhill slopes, they have to encounter about the same amount in
uphill slopes, which are navigated at about 0.7740708353271509 m/s. At that
speed, 53.3 km correspond do about 68899 s or 19 hours. We therefore give
individuals a maximum “sensing” range of 68899 seconds. Individuals will know
about patches within this range from their current patch, but not about patches
further away.

 */
    const MAX_SEASON_RANGE: f64 = 320. / 6. * 1000. / 0.774_070_835_327_150_9;
}
/**

### 4.8 Interaction

> What kinds of interactions among agents are assumed? Are there direct
> interactions in which individuals encounter and affect others, or are
> interactions indirect, e.g., via competition for a mediating resource? If
> the interactions involve communication, how are such communications
> represented? [@grimm2010odd]

 */
pub mod interaction {
    use crate::{emergence, stochasticity, Band, Family, OneYearResources, Parameters};
/**

Agents interact in two different contexts.

THIS USED TO BE:
During the migration phase (part 1 of the step, see [above]), agents avoid
locations where speakers of their language are in the minority, and violently
clash with speakers of other languages. They are attracted by locations where
speakers of their language have large resource stockpiles. These interactions
are described in more detail in [Optimization] and [Adaptation].

*/
    pub fn encounter_maybe_join(family: &mut Family, group: &mut Band, p: &Parameters) -> bool {
        let mut join = true;
        for other_family in &mut group.families {
            if !emergence::will_cooperate(
                &family.culture,
                &other_family.culture,
                p.cooperation_threshold,
            ) {
                join = false;
                let reduction = std::cmp::min(family.effective_size, other_family.effective_size);
                if stochasticity::fight_is_deadly(p.fight_deadliness) {
                    other_family.effective_size -= reduction;
                    if other_family.effective_size == 0 {
                        family.stored_resources += other_family.stored_resources;
                        other_family.stored_resources = OneYearResources::default();
                    }
                }
                if stochasticity::fight_is_deadly(p.fight_deadliness) {
                    family.effective_size -= reduction;
                }
            }
        }
        join
    }
/**

INSTEAD IT IS NOW:
During the migration phase, agents rate their target locations, 

In the resource extraction phase (part 2 of the step, see [above]), agents
compete for limited resources. There is a maximum of resources that can be
extracted from a patch, and the extracted resources are distributed among all
agents in that patch. This distribution is not proportional to the effective
population size: Members of a bigger group of cooperators competing with a
smaller band get an over-proportional part of the total extracted.

*/
    pub fn group_size_effect(raw_group_size: f64) -> f64 {
        raw_group_size.powi(2)
    }
}
/**

*/

/**
### 4.9 Stochasticity

> What processes are modeled by assuming they are random or partly random? Is
> stochasticity used, for example, to reproduce variability in processes for
> which it is unimportant to model the actual causes of the variability? Is it
> used to cause model events or behaviors to occur with a specified frequency?

*/
pub mod stochasticity {
    use crate::Culture;
    /**

    A core random component is the evolution of languages. Our model abstracts away
    from the complexity of real-life language evolution and implements only a
    minimal evolutionary model capable of inheritance and drift. As such, mutations
    happen with a constant rate

    */
    pub fn time_till_next_mutation(culture_mutation_rate: f64) -> u32 {
        fastrand::f64().log(1. - culture_mutation_rate) as u32
    }
    /**
    to random features in ‘linguistic genome’, which is just represented as a binary vector.

    */
    pub fn change_random_bitvec_component(c: &mut Culture) {
        let i = fastrand::usize(0..c.in_memory.len());
        let j = fastrand::u32(0..usize::BITS);
        let flip_j = 1 << j;
        match c.in_memory.get_mut(i) {
            None => {}
            Some(x) => *x ^= flip_j,
        };
    }
    /**
    The linguistic assimilation, where families adopt linguistic features of other
    families they have cooperated with, also relies on stochasticity to choose one
    pivotal family whose linguistic features serve as representative for the whole
    group.

    */
    pub fn select_pivot_family<F>(families: &[F]) -> &F {
        let index = fastrand::usize(0..families.len());
        return families.get(index).unwrap();
    }

    /**
    The other applications of stochasticity serve in the migration process to
    introduce variability compatible with the real world, without having to model
    the underlying causes, of

     - ecological variability, where the availability of resources each season is not deterministic,

       */
    pub fn recover_resources_proportion() -> f64 {
        2.0 * fastrand::f64()
    }
    /**
     - uncertainty in perception when trying to discern whether a potential destination is valuable or not,

    */
    pub fn resource_uncertainty() -> f64 {
        fastrand::f64()
    }
    /**
     - warfare strategies and the deadliness of ambushes, and

    */
    pub fn fight_is_deadly(fight_deadliness: u32) -> bool {
        fastrand::u32(..) < fight_deadliness
    }
    /**
     - who arrives earlier or later at a given migration destination.

    */
    pub fn shuffle<P>(families: &mut Vec<P>) {
        fastrand::shuffle(families);
    }
}

/**
### 4.10 Collectives

> Do the individuals form or belong to aggregations that affect, and are
> affected by, the individuals? How are collectives represented? Is a particular
> collective an emergent property of the individuals, such as a flock of birds
> that assembles as a result of individual behaviors, or is the collective
> simply a definition by the modeler, such as the set of individuals with
> certain properties, defined as a separate kind of entity with its own state
> variables and traits? [@grimm2010odd]

 */
pub mod collectives {
    use crate::{interaction, stochasticity, Band, Family, Parameters};

    /**
    Individual families with compatible languages form ‘bands’ when they are in the
    same patch of land. Families are assigned band membership in random order,
    joining the first band that they encounter which consists only of families with
    a compatible language. As such, membership in bands has a stochastic component:

    With a threshold of 1 and families with languages 00, 01, and 10, there are two
    equally-likely compositons of bands. 01 and 10 are not compatible, because thee
    languages differ in 2 features, and a family with language 00 will join
    whichever of these two groups they encounter first.

    All other properties of a band are derived from the properties of its members,
    but not always in a linear fashion (see eg. [Interaction]). Bands do not persist
    beyond a single time step, they are aggregated after each migration step before
    the resource extraction step to serve as a bookkeeping unit for resource
    extraction and cultural change.

    */
    pub fn cooperate_or_fight<'a>(
        mut families_in_this_location: Vec<&'a mut Family>,
        p: &Parameters,
    ) -> Vec<Band<'a>> {
        stochasticity::shuffle(&mut families_in_this_location);

        let mut bands: Vec<Band<'a>> = vec![];

        for family in families_in_this_location.drain(..) {
            let mut joined_band: Option<&mut Band<'a>> = None;
            for band in &mut bands {
                if interaction::encounter_maybe_join(family, band, p) {
                    joined_band = Some(band);
                    break;
                }
            }
            match joined_band {
                None => {
                    bands.push(Band {
                        culture: family.culture.clone(),
                        families: vec![family],
                    });
                }
                Some(band) => {
                    band.families.push(family);
                }
            }
        }
        bands
            .drain(..)
            .filter_map(|mut band| {
                let pivot_family = stochasticity::select_pivot_family(&band.families);
                band.culture = pivot_family.culture.clone();

                let mut bandsize = 0;
                for family in &mut band.families {
                    bandsize += family.effective_size;
                }
                if bandsize > 0 {
                    Some(band)
                } else {
                    None
                }
            })
            .collect()
    }
}
/**
### 4.11 Observation

> What data are collected from the ABM for testing, understanding, and analyzing
> it, and how and when are they collected? Are all output data freely used, or
> are only certain data sampled and used, to imitate what can be observed in an
> empirical study (‘virtual ecologist’)? [@grimm2010odd]

*/
pub mod observation {
    use crate::{
        Culture, Deserialize, FxHashMap, IntoParallelRefIterator, NodeId, ParallelIterator,
        Parameters, Seasons, Serialize,
    };

    /**
    In order to get an overview over the flow of human migration captured by the
    model and the shape of the language areas, we periodically collect the
    linguistic landscape, i.e. for each language its total population (the sum of
    all families' effective sizes), for each spot.

    */
    pub fn print_population_by_location(
        cultures_by_location: &FxHashMap<NodeId, FxHashMap<Culture, usize>>,
        p: &Parameters,
    ) {
        println!(
            "POPULATION: {:?}",
            cultures_by_location
                .par_iter()
                .map(|(location, cultures)| {
                    let node = &p.dispersal_graph[*location];
                    (
                        node.1,
                        node.2,
                        cultures
                            .iter()
                            .map(|(culture, &count)| (format!("{:x}", culture), count))
                            .collect(),
                    )
                })
                .collect::<Vec<(f64, f64, FxHashMap<String, usize>)>>()
        );
    }

    /**
    Observation parameters govern the logging (and state storing, for inspection and
    resuming) period, in time steps, of the model.

    */
    #[derive(Debug, Serialize, Deserialize, Clone)]
    pub struct Settings {
        pub log_every: Seasons,
        pub log_patch_resources: Seasons,
        pub store_every: Seasons,
        pub statefile: String,
    }
}

/**
## 5. Initialization

> What is the initial state of the model world, i.e., at time $t = 0$ of a
> simulation run? [@grimm2010odd]

> In detail, how many entities of what type are there initially, and what are
> the exact values of their state variables (or how were they set
> stochastically)? [@grimm2010odd]

> Is initialization always the same, or is it allowed to vary among simulations? [@grimm2010odd]

> Are the initial values chosen arbitrarily or based on data? References to
> those data should be provided. [@grimm2010odd]

*/
pub fn very_coarse_dist(x1: f64, y1: f64, x2: f64, y2: f64) -> f64 {
    (x1 - x2).abs() + (y1 - y2).abs()
}

pub fn initialization(p: Parameters, scale: f64) -> State {
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

    let mut patches: FxHashMap<NodeId, Option<Patch>> = FxHashMap::default();
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

    if (p.culture_dimensionality != 1) && (p.culture_dimensionality != 64) {
        panic!(
            "Currently, culture dimensionalities that are not 1 usize/64 bit are not implemented."
        )
    }

    State {
        patches: patches
            .drain()
            .filter_map(|(i, p)| p.map(|q| (i, q)))
            .collect(),
        families: vec![
            Family {
                descendence: String::from("A"),
                location: start1,
                history: vec![],
                seasons_till_next_adult: 4,
                culture: Culture {
                    in_memory: vec![0xffff],
                },

                effective_size: 5,
                number_offspring: 0,
                seasons_till_next_mutation: None,
                stored_resources: OneYearResources::from(10.),
                adaptation: ecology::Ecovector::from(p.minimum_adaptation),
                last_harvest: OneYearResources::default(),
            },
            Family {
                descendence: String::from("F"),
                location: start2,
                history: vec![],
                seasons_till_next_adult: 4,
                culture: Culture { in_memory: vec![0] },

                effective_size: 5,
                number_offspring: 0,
                seasons_till_next_mutation: None,
                stored_resources: OneYearResources::from(10.),
                adaptation: ecology::Ecovector::from(p.minimum_adaptation),
                last_harvest: OneYearResources::default(),
            },
        ],
        t: 0,
        p: p,
    }
}

/**
## 6. Input Data

> Does the model use input from external sources such as data files or other
> models to represent processes that change over time? [@grimm2010odd]

The model currently takes no external data sources into account. As such, it
cannot represent the paleogeographical change in climate, ecology, and geography
(eg. coast lines) that would be necessary to generate a realistic picture of
human migration into the Americas.

*/
mod input {}

/**
## 7. Submodels

> What, in detail, are the submodels that represent the processes listed in
> ‘Process overview and scheduling’? What are the model parameters, their
> dimensions, and reference values? How were submodels designed or chosen, and
> how were they parameterized and then tested? [@grimm2010odd]

*/
pub mod submodels {
    /**
    ### 7.1 Culture

    In order to display evolutionary dynamics, a complex system must at its minimum have
    properties that can be passed on from ancestors to descendants, which undergo variation
    throughout time. For evolution in the classical non-neutral sense, it is furthermore nec-
    essary for the properties to be adaptive, i. e. have an influence on survival and number
    of descendants. The discussion to what extent culture in general is adaptive vs. largely
    neutral is ongoing (). There certainly are cultural domains that show adaptation and have
    effects on survivability, such as climate-appropriate clothing and shelter () or methods of
    subsistence from gathering to agriculture ().

    Of these, agriculture is explicitly outside the scope of the current model, because it
    introduces a feedback loop between culture and carrying capacity which would mask other
    effects. A vast part of the history of humanity went about without agriculture (). But once
    it arose, around 5000 to 10000 years ago independently in different parts of the world (), it
    presumably became a major driver of cultural spreads (farming/language dispersal hypothesis)
    and as such it will be useful to add such effects to a later iteration of this model.
    Certain cultural traits are obviously important for migration processes. Due to lack
    of comparable data on such processes on the family level, we will take an ecology-driven
    probability density to govern knowledge of neighboring patches, and assume this captures
    the direct interactions between the cultural process and the demographic process. We
    describe the motivation of this approximation below in [].

    Beyond traits directly interacting with the ecological niche of a society, there
    is significant literature on the interaction of some very specific cultural
    traits (eg. @Rusch (2014) on altruism and inter-group conflict, Watts et al.
    (2016) on stratified societies and human sacrifice), but very few models that
    consider specific meaniningful cultural traits in agent-based, broad, and
    geographically expansive models. The closest to our goals here may be Hofstede
    et al. (2012). They locate their agents on Hofstede (2001)’s cultural dimensions
    (which are themselves not beyond harsh criticism (McSweeney 2002, Baskerville
    2003)) in order to model negotiations between agents of different cultural
    backgrounds.

    Diverging hypotheses notwithstanding (tone-humidity, elevation-ejectives,
    etc), it seems that a vast number of cultural traits are not directly involved
    in niche adaptation. This does not mean individual traits do not exhibit
    non-neutral evolutionary dynamics, as shown eg. by Kandler & Shennan (2013). But
    a neutral model for cultural change appears appropriate nonetheless.

    > FIXME: Kandler and others have also argued that mixing a large amount of
    > non-neutral evolutionary processes is not well approximated by generally
    > neutral evolution, so this needs a slightly stronger case.

    Neutral abstract models for culture have in the past been considered for various
    purposes (Komarova et al. 2001). In such models, culture tends to be modeled as a binary
    vector (Fogarty et al. 2017, Pascual et al. 2020). The number of dimension M of this culture
    vector range between M = 6 (Pascual et al. 2020) and M = 10 [@delcastillo2013modeling],
    with empirical reasons cited for ??? [@?] and theoretical reasons for ??? [@?] or M > N for the (effec-
    tive) population size N (Fogarty et al. 2017). Vectors are compared using the Hamming
    distance (after, measuring the number of mismatches between the two vectors) in most
    cases (Fogarty et al. 2017, Pascual et al. 2020).

    Other options exist. For instance, [@Barceló et al. (2014); @Barceló et al.
    (2015)] use integer vectors with values between 1 and 6 reflecting importance,
    and ‘a multidimensional weighted Euclidean distance based on the extension of
    the cosine similarity measure for vectors’ (2014), with weights that ‘roughly
    imitate[] the results of a factorial analysis of individual beliefs’ (2014). It
    has not been shown that this approach to modeling culture improves realism or
    interpretability, so we take a neutral evolution model with binary vectors and
    an unweighted Hamming distance.

    */
    pub mod culture {
        use crate::stochasticity;
        use crate::Culture;

        /** At a given low rate, mutate a random feature of a language.

        FIXME: Conceptually, it would make sense to store the
        seasons_till_next_mutation inside the Culture struct, but then we need
        to be careful to not copy it when we clone those objects.
         */
        pub fn mutate(
            family_seasons_till_next_mutation: &mut Option<u32>,
            family_culture: &mut Culture,
            target_culture: &Culture,
            culture_dimensionality: usize,
            culture_mutation_rate: f64,
        ) {
            *family_culture = target_culture.clone();
            match family_seasons_till_next_mutation {
                None => {
                    *family_seasons_till_next_mutation = Some(
                        stochasticity::time_till_next_mutation(culture_mutation_rate),
                    );
                    mutate(
                        family_seasons_till_next_mutation,
                        family_culture,
                        target_culture,
                        culture_dimensionality,
                        culture_mutation_rate,
                    );
                }
                Some(0) => {
                    stochasticity::change_random_bitvec_component(family_culture);
                    *family_seasons_till_next_mutation = Some(
                        stochasticity::time_till_next_mutation(culture_mutation_rate),
                    );
                }
                Some(k) => {
                    *family_seasons_till_next_mutation = Some(*k - 1);
                }
            }
        }
    }

    pub mod family_lifecycle {
        use crate::{adaptation, sensing};
        use crate::{Culture, Family, NodeId, OneYearResources, Parameters, Patch};
        use dashmap::DashMap;
        use rayon::prelude::*;
        use rustc_hash::{FxHashMap, FxHasher};

        pub fn procreate_and_migrate(
            families: &mut Vec<Family>,
            storage: &DashMap<
                NodeId,
                FxHashMap<Culture, (usize, OneYearResources)>,
                std::hash::BuildHasherDefault<FxHasher>,
            >,
            p: &Parameters,
            patches: &mut DashMap<NodeId, Patch>,
            cc: &DashMap<NodeId, usize, std::hash::BuildHasherDefault<FxHasher>>,
            nearby_cache: &mut DashMap<
                NodeId,
                Vec<NodeId>,
                std::hash::BuildHasherDefault<FxHasher>,
            >,
        ) -> Vec<Family> {
            let cache = std::sync::Arc::new(DashMap::default());

            families
                .par_iter_mut()
                .filter_map(|mut family| {
                    let nearby = nearby_cache
                        .entry(family.location)
                        .or_insert_with(|| sensing::nearby_locations(family.location, p));
                    // println!("{:?}", *nearby);

                    let (destination, cost) = adaptation::decide_on_moving(
                        family, storage, &nearby, patches, &cache, p, cc,
                    );
                    // println!("Family {:} moved to {:?}", family.descendence, destination);

                    family.history.push((
                        family.location,
                        family.last_harvest / family.effective_size as f64,
                    ));
                    family.stored_resources += family.last_harvest;
                    family.last_harvest = OneYearResources::default();
                    family.location = destination;
                    family.stored_resources -= cost;

                    match maybe_procreate(family, p.season_length_in_years) {
                        None => None,
                        // In terms of scheduling a new family can (and if possible
                        // should) move immediately when created. This behaviour is
                        // taken from del Castillo (2013).
                        Some(mut descendant) => {
                            // println!("Family {:} split off family {:}", descendant.descendence, family.descendence);
                            let (destination, cost) = adaptation::decide_on_moving(
                                &descendant,
                                storage,
                                // TODO: Check that the center is at index 0
                                &nearby[1..],
                                patches,
                                &cache,
                                p,
                                cc,
                            );
                            descendant.location = destination;
                            descendant.stored_resources -= cost;
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
                let resources_they_take =
                    family.stored_resources * 2.0 / family.effective_size as f64;
                family.stored_resources -= resources_they_take;
                family.effective_size -= 2;
                let time_to_grow_adult = (15. / season_length_in_years) as u32 + 1;
                Some(Family {
                    descendence: format!("{}:{:}", family.descendence, family.number_offspring),
                    location: family.location,
                    history: family.history[family.history.len()
                        - std::cmp::min(family.history.len(), time_to_grow_adult as usize)
                        ..family.history.len()]
                        .into(),
                    seasons_till_next_adult: time_to_grow_adult,
                    culture: family.culture.clone(),

                    effective_size: 2,
                    number_offspring: 0,
                    seasons_till_next_mutation: None,
                    stored_resources: resources_they_take,
                    adaptation: family.adaptation,
                    last_harvest: OneYearResources::default(),
                })
            }
        }

        /**
        According to [@goodman1985menarche], the spacing of births where infant survives
        until birth of next sibling is $2.85±1.35$. It has been hypothesized [@?] that
        the interval may be lower in populations at the pop cap, so we take μ-σ.
        According to [@volk2013infant], the mean child mortality rate
        (cumulative probability of dying prior to approximate sexual maturity at age 15)
        of modern hunter-gatherers is 48.8%, so the mean interval of new adults is the
        mean interval of children, divided by one minus that rate.

        */
        const YEARS_BETWEEN_ADULTS: f64 = (2.85 - 1.35) / (1. - 0.488);

        pub fn maybe_grow(family: &mut Family, season_length_in_years: f64) {
            if family.seasons_till_next_adult == 0 {
                family.effective_size += 1;
                family.seasons_till_next_adult =
                    (YEARS_BETWEEN_ADULTS / season_length_in_years) as u32;
            } else {
                family.seasons_till_next_adult -= 1;
            }
        }

        // Find at least a proxy for children using resources, or model children explicitly
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
        use crate::{interaction, stochasticity, OneYearResources, Parameters};

        /**

                                                                                                                                                                        This function models how much a family contributes to the group effort of
                                                                                                                                                                        harvesting a particular ecoregion. The family's contribution is given by its
                                                                                                                                                                        size, weighted with how well they know the ecoregion.

                                                                                                                                                                        */
        /**
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
            // TODO: This means if there are multiple ecoregions, we forget more about the ones we are not in.
            family.adaptation = family.adaptation * 0.95;
            family.adaptation.entries[ecoregion] += 0.05;
            let contribution = family.effective_size as f64 * family.adaptation[ecoregion];
            *sum_effort += contribution;
            (&mut family.last_harvest, contribution)
        }

        pub fn adjust_culture(group: &mut crate::Band, p: &Parameters) {
            for family in &mut group.families {
                crate::submodels::culture::mutate(
                    &mut family.seasons_till_next_mutation,
                    &mut family.culture,
                    &group.culture,
                    p.culture_dimensionality,
                    p.culture_mutation_rate,
                );
            }
        }

        pub fn group_contribution<'a>(
            group: &'a mut crate::Band,
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
        pub fn exploit_patch<P>(mut groups: Vec<crate::Band>, mut patch: P, p: &Parameters)
        where
            P: core::ops::DerefMut<Target = crate::Patch>,
        {
            for (ecoregion, (res, max_res)) in &mut patch.resources {
                let mut total_squared_groups_effort = 1.0;
                let mut work = groups
                    .iter_mut()
                    .map(|group| {
                        let mut sum_effort = 0.0;
                        let contributions: Vec<_> =
                            group_contribution(group, *ecoregion, &mut sum_effort);
                        total_squared_groups_effort += interaction::group_size_effect(sum_effort);
                        (sum_effort, contributions)
                    })
                    .collect::<Vec<(f64, Vec<(&mut OneYearResources, f64)>)>>();
                // println!("{:?} from {:?}", total_squared_groups_effort, work);
                work.drain(..)
                    .for_each(|(group_effort, mut contributions)| {
                        assert!(group_effort > 0.);
                        // println!("group effort: {:?}", group_effort);
                        let group_harvest = std::cmp::min(
                            *res * interaction::group_size_effect(group_effort)
                                / total_squared_groups_effort,
                            p.maximum_resources_one_adult_can_harvest
                                * p.season_length_in_years
                                * group_effort,
                        );
                        assert!(group_harvest > OneYearResources::from(0.));
                        // println!("group harvest: {:?}", group_harvest);
                        let coefficient = std::cmp::min(
                            group_harvest / group_effort,
                            OneYearResources::from(2.5),
                        );
                        assert!(coefficient > OneYearResources::from(0.));
                        // println!("coefficient: {:?}", coefficient);
                        for (harvest, contribution) in &mut contributions {
                            **harvest = coefficient * p.season_length_in_years * *contribution;
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
                    patch_max_resources
                        * resource_recovery_per_season
                        * stochasticity::recover_resources_proportion()
                        + *patch_resources,
                )
            }
        }
    }

    pub mod parameters {
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
    }
}

pub fn store_state(state: State, statefile: String) -> Result<(), String> {
    let file = match File::create(statefile) {
        Ok(f) => f,
        Err(_) => return Err("Could not create state file".to_string()),
    };
    match bincode::serialize_into(file, &state) {
        Ok(_) => Ok(()),
        Err(_) => Err("Failed to store state.".to_string()),
    }
}

pub fn run(mut s: State, max_t: Seasons, o: &observation::Settings) {
    let mut stored_resources: DashMap<_, _, std::hash::BuildHasherDefault<FxHasher>> =
        DashMap::default();
    let mut cache: DashMap<_, _, std::hash::BuildHasherDefault<FxHasher>> = DashMap::default();
    loop {
        stored_resources = step(
            &mut s.families,
            &mut s.patches,
            &stored_resources,
            &s.p,
            s.t,
            &mut cache,
            o,
        );

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
