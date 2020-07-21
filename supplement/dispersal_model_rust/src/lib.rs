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
split, driven only by the necessary interactions between humans.

The summary statistics of this phylogeny (in particular diversification rates)
are to be compared to values known from language evolution. The model is
structured to be easily applied to study the history of the settlement of the
Americas at a later time. It would require paleoclimate data and interpretations
of past ecoregions to produce results that can be compared to that history.
*/

// TODO: Dear Rustacean, I know that my use of documentation comments is
// hazardous, because they land in the generated documentation in a different
// order, or attached to things that are only accidentally coming afterwards.

// Load useful modules

use rand::prelude::*;
use rayon::prelude::*;
use std::collections::{HashMap,BTreeSet};
use rand_distr::StandardNormal;

mod debug;
pub mod ecology;

use submodels::parameters::Parameters;

/**

# 2. Entities, state variables, and scales

The model consists of agents interacting on a hexagonal discrete global grid in
discrete time. One time step is supposed to model a season with a duration of
half a year.
 */

pub type HalfYears = u32;

/**
Whereever possible, resources are measured in kcal (in SI units: 1 kcal = 4.184 kJ)
 */

pub type KCal = f32;

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
pub type NodeId = usize;

/*
Each contains one or more ecoregions, according to the geographical distribution
of ecoregions in the real world. For each ecoregion, the patch has a maximum and
current availability of resources, measured in kcal available over the course of
half a year. These numbers are estimated from population densities following
[@tallavaara2018resource] in conjunction with the ecoregion-covered hex areas in
km².

(In the database containing the pre-processed real-world data, the hex area is
given in cosine-rescaled 15" squares, and multiplicative constants are used to
translate those numbers into km².)
 */


/**
In the current state of the model, the resource maximum is constant throughout
the simulation, but future work might add seasonality (with every other time
step corresponding to the more plentiful half of the year and the others to the
scarcer half) or random variation (eg. estimated from the observed random
variation of resources within that ecoregion). By nature, the simulation invites
the extension to include paleoclimate data, but that only becomes relevant once
it shows fundamentally reasonable dynamics.
 */

pub struct Patch {
    /// For every local ecoregion, the tuple of resources currently accessible
    /// in this patch, and the maximum possible resources available in a time
    /// step, both in kcal / season.
    resources: HashMap<usize, (KCal, KCal)>,
}

/**
## 2.2 Families

The main decision-making agents of the simulation are families. Families can
migrate between cells and form links to other families in the context of
cooperation to extract resources.
 */

pub struct Family {
    /// The agent's history of decendence, also serving as unique ID.
    descendence: String,
    /// The effective size of the family in number of adults. One adult is
    /// assumed to consume the same amount of food/energy and to contribute the
    /// same labor to foraging as any other adult.
    effective_size: usize,
    /// The family's shared culture.
    culture: Culture,
    /// The current location of the agent.
    location: NodeId,
    /// The previous locations of the agent. This is useful for some bits of
    /// analysis, but agents are also guaranteed to have knowledge about their
    /// recent locations, while knowlede beyond their own history depends on
    /// random factors. For their historical locations, the family also knows
    /// the payoff they gained there, so we keep these two bits of information
    /// together in this place.
    location_history: Vec<(NodeId, KCal)>,
    /// The amount of stored resources, in kcal, the family has access to
    /// without going foraging.
    stored_resources: KCal,
    /// Adaptation to local ecoregions is represented as a vector with values
    /// between 0.0 (completely unknown) and 1.0 (perfectly familiar) for any
    /// ecoregion the family might encounter.
    adaptation: ecology::Ecovector,
    /// The number of seasons to wait until the next child (reset to 2 when
    /// starvation happens)
    seasons_till_next_child: HalfYears,

    /// For bookkeeping purposes (eg. generating descendant's ‘descendence’
    /// values), keep track of the number of offspring families this family has
    /// spawned so far.
    number_offspring: u16,
    /// A bookkeeping quantity. Instead of mutation happening in each time step
    /// with a tiny probability, the same distribution of times between
    /// mutations is generated by drawing a number of season from a geometric
    /// distribution and counting it down by one each time step, which is useful
    /// because random number generation is computationally expensive.
    seasons_till_next_mutation: Option<HalfYears>,
}

impl PartialEq for Family {
    fn eq(&self, other: &Self) -> bool {
        self.descendence == other.descendence
    }
}

/**
Families in the same location with compatible cultures can cooperate to improve
their success at extracting resources. (Following XXX, cooperation could also
mean sharing resources between different families. To simplify the model, this
implementation does not contain that effect.) The cooperative groups formed by
cooperating families are higher-level agents created ad-hoc in each time step.
They do not persist or have effect beyond a single time step. Resource
exploitation happens at the level of the cooperative group and is distributed to
the individual families after the fact.
 */

/**
## 2.3 Cultures

Every family has a culture. These are very abstract and vastly simplified, due
to the lack of quantitative data on cultural evolution in a framework comparable
to the one used for this study. Based on the need to have a data structure that
supports random drift, a low but positive chance of back-mutation, and a large
number of equally similar cultures, we describe culture using a binary vector.
Computers make natural use of the binary representation integer numbers, so a
binary vector is equivalent to an unsigned integer of sufficient size, which is
faster to use in computations and more efficient to store.
 */
type Culture = u64;

/**
## 2.4 State

The associations between gridcells and patches and between gridcells and the
families located there (stored in the Family's `location`) are the core of the
model state. The state also tracks the time, measured in time steps
corresponding to half a year each, since the start of the simulation, and stores
a copy of the model parameters.
 */
pub struct State {
    /// The patches of the model, indexed by their address according to the H3
    /// discrete global grid system. This set is fixed, although the properties
    /// of the patches may change with time.
    patches: HashMap<NodeId, Patch>,
    /// The agents (families) currently active in the simulation. This vector
    /// changes over time as new families are added and old families die out.
    families: Vec<Family>,
    /// The current time step, in half years since start of the simulation.
    t: HalfYears,
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
    patches: &mut HashMap<NodeId, Patch>,
    knowledge: HashMap<NodeId, Knowledge>,
    p: &Parameters,
    t: HalfYears,
) -> (HashMap<NodeId, Vec<Family>>, HashMap<NodeId, Knowledge>)
{
    let mut families_by_location = step_part_1(families, patches, knowledge, p, t);
    let new_knowledge = step_part_2(&mut families_by_location, patches, p);
    (families_by_location, new_knowledge)
}

/**
The first part focuses on the individual families, which shrink, grow, die,
split, and move. It constructs the mapping of families at the end of the season,
grouped by their location after potential moves. Because the movement of a
family depends on the distribution of othe families at he start of the season,
it can happen entirely in parallel.
 */

fn step_part_1(
    families: &mut Vec<Family>,
    patches: &mut HashMap<NodeId, Patch>,
    knowledge: HashMap<NodeId, Knowledge>,
    p: &Parameters,
    t: HalfYears,
) -> HashMap<NodeId, Vec<Family>>
{
    let mut families_by_location: HashMap<NodeId, Vec<Family>> = HashMap::new();
    {
        // For reporting only
        let mut cultures_by_location: HashMap<NodeId, HashMap<Culture, usize>> = HashMap::new();

        for family in families.iter_mut() {
            family.adaptation = family.adaptation * 0.95;
            if submodels::family_lifecycle::use_resources_and_maybe_shrink(
                &mut family.effective_size,
                &mut family.stored_resources,
                p,
            ) {
                family.seasons_till_next_child = std::cmp::max(family.seasons_till_next_child, 2)
            }

            submodels::family_lifecycle::maybe_grow(family);

            let cultures = cultures_by_location
                .entry(family.location)
                .or_insert_with(HashMap::new);
            let counter = cultures.entry(family.culture).or_insert(0);
            *counter += family.effective_size;
        }

        if t % 20 == 0 {
            let l_c = cultures_by_location.clone();
            observation::print_gd_cd(l_c, p);
            let l_c = cultures_by_location;
            observation::print_population_by_location(l_c, p);
        }
    }

    let mut rng = rand::thread_rng();
    families.shuffle(&mut rng);

    for mut family in families.drain(..) {
        if submodels::family_lifecycle::is_moribund(&family) {
            family.effective_size = 0;
            // Stop tracking this family, they are dead.
            continue;
        }

        let nearby = sensing::nearby_locations(family.location, p);
        println!("{:?}", nearby);

        match submodels::family_lifecycle::maybe_procreate(&mut family) {
            None => {}
            // In terms of scheduling a new family can (and if possible
            // should) move immediately when created. This behaviour is
            // taken from del Castillo (2013). Also, a descendant family
            // will move directly before their progenitor.
            Some(descendant) => {
                let d = adaptation::decide_on_moving(
                    &descendant,
                    nearby.iter().filter_map(
                        |(i, d)| if *i == family.location {
                            None
                        } else {
                            Some((*i, *d, patches.get(i)?, knowledge.get(i)))
                        }),
                    true, p);
                let destination = d.unwrap_or(family.location);

                families_by_location
                    .entry(destination)
                    .or_insert(vec![])
                    .push(descendant);
            }
        }
        let d = adaptation::decide_on_moving(
            &family,
            nearby.iter().filter_map(
                |(i, d)|
                Some((*i, *d, patches.get(i)?, knowledge.get(i)))
            ),
            false, p);

        // Update cultures_by_location
        let destination = d.unwrap_or(family.location);
        families_by_location
            .entry(destination)
            .or_insert(vec![])
            .push(family);
    }
    families_by_location
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
    families_by_location: &mut HashMap<NodeId, Vec<Family>>,
    patches: &mut HashMap<NodeId, Patch>,
    p: &Parameters) -> HashMap<NodeId, Knowledge>
{
    let mut rng = rand::thread_rng();
    let mut knowledge = HashMap::new();

    for (patch_id, families) in families_by_location {
        let mut sum_normalized_payoff = 0.0;
        let mut cultures: BTreeSet<Culture> = BTreeSet::new();

        assert!(!families.is_empty());

        let patch = match patches.get_mut(patch_id) {
            None => { continue },
            Some(p) => { p }
        };

        let mut groups = collectives::cooperatives(
            families.iter_mut().collect(), p);
        let n_groups = groups.len() as KCal;

        let mut unextracted = 0.0;
        for (i, (res, res_max)) in patch.resources.iter_mut() {
            let mut sum_extracted = 0.0;
            let actual_contributions: Vec<(_, f32, f32)> = groups.iter_mut().map(|group| {
                let mut raw_contribution = 0.0;
                let target_culture = group.culture;
                cultures.insert(target_culture);
                let families: Vec<(&mut f32, f32)> = group.families.iter_mut().map(|family| {
                    family.adaptation.entries[*i] += 0.05;
                    submodels::culture::mutate_culture(
                        &mut family.seasons_till_next_mutation,
                        &mut family.culture, target_culture,
                        p.culture_dimensionality, p.culture_mutation_rate);
                    let contribution = family.effective_size as f32 * family.adaptation[*i];
                    raw_contribution += contribution;
                    (&mut family.stored_resources, contribution)
                }).collect();
                let actual_contribution: f32 =
                    interaction::effective_labor_through_cooperation(
                        raw_contribution +
                            (p.payoff_std *
                             rng.sample::<f32, _>(StandardNormal) *
                             raw_contribution.powf(0.5)),
                        p.cooperation_gain);
                assert!(actual_contribution.is_normal());
                sum_extracted += actual_contribution;
                (families, raw_contribution, actual_contribution)
            }).collect();
            let actual_payout = f32::min(*res, sum_extracted);
            unextracted += *res - actual_payout;
            *res -= actual_payout;
            for (families, raw_contribution, actual_contribution) in actual_contributions {
                let actual_extracted = actual_payout * actual_contribution / sum_extracted;
                sum_normalized_payoff += actual_extracted / raw_contribution;
                for (family_res, contribution) in families {
                    let actual_returns = actual_extracted * contribution / raw_contribution;
                    *family_res += actual_returns;
                }
            }
            submodels::ecology::recover(
                res, res_max, p.resource_recovery, p.accessible_resources);
        }
        knowledge.insert(*patch_id, Knowledge {
            local_cultures: cultures,
            normalized_payoff: sum_normalized_payoff / n_groups,
            leftover: unextracted,
        });
    }
    knowledge
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
            let l1h3 = p.dispersal_graph[*l1].0;
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
                                .or_insert(0)
                                += count1 * (count2 + 1) / 2;
                            break;
                        }
                        let cd = submodels::culture::distance(*c1, *c2);
                        *gd_cd
                            .entry(gd)
                            .or_insert_with(HashMap::new)
                            .entry(cd)
                            .or_insert(0)
                            += count1 * count2;
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
     */

    /**
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

Agents also adapt to local environment: TODO
*/

pub struct Knowledge {
    local_cultures: BTreeSet<Culture>,
    normalized_payoff: KCal,
    leftover: KCal,
}

const WALKING_TIME_PER_SEASON: f64 = 365.24219 * 0.5 * NORM;
const COST_PER_WALKING_TIME: f64 =
    // Triple: opportunity cost from not foraging, and cost from doing likely
    // heavy labor instead.
    3. / WALKING_TIME_PER_SEASON;

mod adaptation {
    use crate::*;

    pub fn decide_on_moving<'a, KD>(
        family: &'a Family,
        kd: KD,
        avoid_stay: bool,
        p: &Parameters,
    ) -> Option<NodeId>
    where
        KD: Iterator<Item = (NodeId, f64, &'a Patch, Option<&'a Knowledge>)>,
    {
        let threshold = if avoid_stay {
            0.
        } else {
            p.time_step_energy_use * p.evidence_needed +
                match family.location_history.last() {
                    None => 0.0,
                    Some((_, k)) => *k
                }
        };

        objectives::best_location(
            family.location_history.iter().map(|(i, j)| (*i, *j)).chain(kd.map(
                |(i, d, l, k)| {
                    let movement_cost = (d * COST_PER_WALKING_TIME) as f32 * p.time_step_energy_use ;
                    let mut tot_res = 0.;
                    let mut family_res = 0.;
                    let sizef = family.effective_size as f32;
                    for (j, (res, _)) in &l.resources {
                        tot_res += res;
                        family_res += res * family.adaptation[*j];
                    };
                    (i,
                     family_res / tot_res * match k {
                         None =>
                             f32::min(interaction::effective_labor_through_cooperation(sizef, p.cooperation_gain) / sizef,
                                      l.resources.values().map(|(res, _)| res).sum()),
                         Some(k) => if k.local_cultures.iter().any(
                             |c| emergence::similar_culture(
                                 *c, family.culture, p.cooperation_threshold)) {
                             f32::min(k.normalized_payoff, k.leftover)
                         } else {
                             f32::min(interaction::effective_labor_through_cooperation(sizef, p.cooperation_gain) / sizef,
                                      k.leftover)
                         }} - movement_cost)})),
            threshold
        ).or(Some(family.location))
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
    pub fn best_location<Pair, N, I>(
        kd: Pair,
        threshold: N,
    ) -> Option<I>
    where
        Pair: Iterator<Item = (I, N)>,
        N: PartialOrd + Default + std::fmt::Debug,
        I: Default + std::fmt::Display,
    {
        let mut rng = rand::thread_rng();

        // This variable `n_best_before` is used to randomly draw between
        // several equally-optimal options. It counts the number of best options
        // encountered so far.
        let mut n_best_before = 0;
        let mut target: Option<I> = None;
        let mut max_gain = N::default();

        for (location, expected_gain) in kd {
            if expected_gain >= max_gain {
                if expected_gain < threshold {
                    continue;
                }
                if expected_gain == max_gain {
                    n_best_before += 1;
                    if rng.gen_range(0, n_best_before + 1) < n_best_before {
                        continue;
                    }
                } else {
                    n_best_before = 0
                }
                target = Some(location);
                max_gain = expected_gain;
            }
        }
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

mod learning {
}

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
mod prediction {
    

}

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

 */

const NORM: f64 = 8. * 60. * 60.; // seconds

mod sensing {
    use crate::*;

    /// Individuals know about nearby locations. The exploration is not
    /// explicitly modelled.

    pub fn nearby_locations(
        location: NodeId,
        p: &Parameters,
    ) -> Vec<(NodeId, f64)> {
        let mut rng = rand::thread_rng();
        movementgraph::bounded_dijkstra(
            &p.dispersal_graph,
            location,
            NORM * 12., //4 complete days, or 12 days of travel
            |e| *petgraph::csr::EdgeReference::weight(&e)
        )
            .iter()
            .filter_map(
                |(n, v)| {
                println!("{:}, {:}", n, v);
                if rng.gen::<f64>() < 1./(2. + v / NORM){
                    Some((*n, *v))
                } else {
                    None
                }})
            .collect()
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

    pub fn effective_labor_through_cooperation(n_cooperators: f32, cooperation_gain: f32) -> f32 {
        // 1. + (n_cooperators - 1.).powf(cooperation_gain) / 10.
        n_cooperators.powf(1. + cooperation_gain)
    }
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
mod collectives {
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
        let mut groups: Vec<Cooperative> = vec![];

        for family in families_in_this_location.drain(..) {
            let mut joined_group: Option<&mut Cooperative> = None;
            for group in groups.iter_mut() {
                let mut join = true;
                for other_family in group.families.iter() {
                    if !emergence::similar_culture(family.culture, other_family.culture, p.cooperation_threshold) {
                        join = false;
                        break;
                    }
                }
                if join {
                    joined_group = Some(group);
                    break;
                }
            }
            let _size = family.effective_size as f32;
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
mod observation {
    use crate::*;

    /**
    Do get an overview over the flow of human migration captured by the model,
    we collect the population count, i.e. the sum of all families' effective
    sizes, for each spot in each time step.
    */
    pub fn print_population_by_location(
        cultures_by_location: HashMap<NodeId, HashMap<Culture, usize>>,
        p: &Parameters,
    ) {
        println!(
            "POPULATION: {:?}",
            cultures_by_location
                .par_iter()
                .map(|(k, v)| {
                    let g = hexgrid::geo_coordinates(p.dispersal_graph[*k].0);
                    (g.lon, g.lat, v.iter().map(|(_, c)| c).sum())
                })
                .collect::<Vec<(f64, f64, usize)>>()
        );
    }

    /**
    The major focus of the model, however, is on the cultural areas and
    networks. To asses the main outcome of the model, we collect the geographic
    vs. cultural distances for every pair of individuals within 5 years of each other. Due to the
    computational effort and high autocorrelation of the outcomes, we collect
    this data only about once per generation, i.e. every 60 time steps.
    */
    pub fn print_gd_cd(
        cultures_by_location: HashMap<NodeId, HashMap<Culture, usize>>,
        p: &Parameters) {
        let gd_cd: HashMap<i32, HashMap<u32, usize>> =
            emergence::cultural_distance_by_geographical_distance(&cultures_by_location, 5 * 5 * 2, p);
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

pub fn initialization(
    p: &Parameters) -> Option<State>
{
    let graph = &p.dispersal_graph;

    let mut start1: NodeId = 0;
    let mut start2: NodeId = 0;
    let mut start1d = very_coarse_dist(graph[0].1, graph[0].2, -159.873, 65.613);
    let mut start2d = very_coarse_dist(graph[0].1, graph[0].2, -158.2718, 60.8071);

    for i in 0..graph.node_count() {
        let longitude = graph[i].1;
        let latitude = graph[i].2;
        let d1 = very_coarse_dist(longitude, latitude, -159.873, 65.613);
        let d2 = very_coarse_dist(longitude, latitude, -158.2718, 60.8071);
        if d1 < start1d {
            start1 = i;
            start1d = d1;
        }
        if d2 < start2d {
            start2 = i;
            start2d = d2;
        }
    }
    let _area: f32 = libh3::hex_area_km_2(5) as f32;

    println!("Starts: {:}, {:}", start1, start2);

    let mut patches: HashMap<NodeId, Option<Patch>> = HashMap::new();
    let mut new_patches: std::collections::BinaryHeap<NodeId> =
        std::collections::BinaryHeap::new();
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

        println!("Resources at {:} ({:}, {:}): {:?}",
                 next, longitude, latitude, ecoregions);
        patches.insert(
            next,
            Some(Patch { resources: ecoregions.iter().map(
                |(&i, &j)| (i, (j * p.time_step_energy_use, j * p.time_step_energy_use))).collect() }),
        );
        for q in graph.neighbors_slice(next) {
            if patches.contains_key(q) {
                continue;
            }
            new_patches.push(*q);
        }
    }
    let f1 = Family {
        descendence: String::from("A"),
        location: start1,
        location_history: vec![],
        seasons_till_next_child: 4,
        culture: 0b000_000_000_000_000,

        effective_size: 5,
        number_offspring: 0,
        seasons_till_next_mutation: None,
        stored_resources: p.time_step_energy_use * 10.,
        adaptation: ecology::Ecovector::default(),
    };
    let f2 = Family {
        descendence: String::from("F"),
        location: start2,
        location_history: vec![],
        seasons_till_next_child: 4,
        culture: 0b111_111_111_111_111,

        effective_size: 5,
        number_offspring: 0,
        seasons_till_next_mutation: None,
        stored_resources: p.time_step_energy_use * 10.,
        adaptation: ecology::Ecovector::default(),
    };

    Some(State {
        patches: patches
            .drain()
            .filter_map(|(i, p)| match p {
                None => None,
                Some(q) => Some((i, q)),
            })
            .collect(),
        families: vec![
            f1
            // f2,
        ],
        t: 0,
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
            (c1 ^ c2).count_ones()
        }

        pub fn mutate_culture(
            family_seasons_till_next_mutation: &mut Option<u32>,
            family_culture: &mut Culture, target_culture: Culture,
            culture_dimensionality: u8, culture_mutation_rate: f64)
        {
            *family_culture = target_culture;
            match family_seasons_till_next_mutation {
                None => {
                    *family_seasons_till_next_mutation =
                        Some(random::<f64>().log(1. - culture_mutation_rate) as u32);
                    mutate_culture(family_seasons_till_next_mutation,
                                   family_culture, target_culture,
                                   culture_dimensionality, culture_mutation_rate);
                },
                Some(0) => {
                    let i: u8 = rand::thread_rng().gen_range(0, culture_dimensionality);
                    *family_culture ^= 1 << i;
                    *family_seasons_till_next_mutation =
                        Some(random::<f64>().log(1. - culture_mutation_rate) as u32);
                },
                Some(k) => {
                    *family_seasons_till_next_mutation = Some(*k - 1);
                }
            }
        }

    }

    pub mod family_lifecycle {
        use crate::{Family, KCal, Parameters};

        pub fn resources_at_season_end(resources: KCal, size: usize, p: &Parameters) -> KCal
        {
            let mut resources_after: KCal = resources - (size as f32) * p.time_step_energy_use;
            if resources_after > 0. {
                resources_after *= 1. - p.storage_loss;
            }
            resources_after
        }

        pub fn maybe_procreate(family: &mut Family) -> Option<Family> {
            if family.effective_size < 10 {
                None
            } else {
                family.number_offspring += 1;
                family.effective_size -= 2;
                Some(Family {
                    descendence: format!("{}:{:}", family.descendence, family.number_offspring),
                    location: family.location,
                    location_history: vec![],
                    seasons_till_next_child: 12 * 2,
                    culture: family.culture,

                    effective_size: 2,
                    number_offspring: 0,
                    seasons_till_next_mutation: None,
                    stored_resources: 0.,
                    adaptation: family.adaptation,
                })
            }
        }

        pub fn maybe_grow(family: &mut Family) {
            if family.seasons_till_next_child == 0 {
                family.effective_size += 1;
                family.seasons_till_next_child = 2;
            } else {
                family.seasons_till_next_child -= 1;
            }
        }
        pub fn use_resources_and_maybe_shrink(
            size: &mut usize,
            resources: &mut KCal,
            p: &Parameters,
        ) -> bool
        {
            let mut has_shrunk = false;
            while resources_at_season_end(*resources, *size, p) < 0. && *size > 0 {
                *size -= 1;
                has_shrunk = true;
            }
            *resources = resources_at_season_end(*resources, *size, p);
            has_shrunk
        }

        pub fn is_moribund(family: &Family) -> bool {
            family.effective_size < 2
        }
    }

    pub mod ecology {
        use crate::KCal;

        /**
        Recover a patch, mutating its resources

        The recovery assumes exponential (or geometric, it's a discrete-time
        model) growth of resources in a patch up to a given maximum.

        The patch's `max_resouces` represents the maximum available for human
        extraction. This is a proportion `accessible_resources` of the total
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
        recover(&mut resources, &2.0, 0.5, 0.25);
        assert_eq!(resources, 1.4375);
        ```

        */
        pub fn recover(
            patch_resources: &mut KCal,
            patch_max_resources: &KCal,
            resource_recovery: f32,
            accessible_resources: f32)
        {
            if *patch_resources < *patch_max_resources {
                let inaccessible = patch_max_resources * (1. - accessible_resources) / accessible_resources;
                let underlying_resources = *patch_resources + inaccessible;
                *patch_resources += underlying_resources
                    * resource_recovery
                    * (1. - underlying_resources * accessible_resources/ patch_max_resources);
            }
            assert!(patch_resources.is_normal(),
                    "The recovery of {} from {} towards {} ({}% of total) was not finite",
                    resource_recovery,
                    *patch_resources,
                    *patch_max_resources,
                    accessible_resources * 100.
            );
        }
    }

    pub mod parameters {
        use crate::KCal;
        use crate::movementgraph::MovementGraph;

        pub struct Parameters {
            pub attention_probability: f32,
            pub time_step_energy_use: KCal,
            pub storage_loss: f32,
            pub resource_recovery: f32,
            pub culture_mutation_rate: f64,
            pub culture_dimensionality: u8,
            pub cooperation_threshold: u32,
            pub cooperation_gain: f32,
            pub accessible_resources: f32,
            pub evidence_needed: f32,
            pub payoff_std: f32,

            pub dispersal_graph: MovementGraph,
        }
    }
}


pub fn run(mut s: State, p: Parameters, max_t: HalfYears) {
    let mut knowledge = HashMap::new();
    loop {
        let (families_by_location, k) = step(&mut s.families, &mut s.patches, knowledge, &p, s.t);
        knowledge = k;
        for (location, families) in families_by_location {
            for mut family in families {
                family.location_history.push((family.location, family.stored_resources));
                family.location = location;
                s.families.push(family);
            }
        }
        if s.families.is_empty() {
            println!("Died out");
            break;
        }
        s.t += 1;
        if s.t > max_t {
            println!("Ended");
            break
        }
    }
}
