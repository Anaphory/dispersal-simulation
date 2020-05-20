/*!
Model Description {#odd}
========================

This model description follows the ODD (Overview, Design concept, Details)
protocol [@grimm2006standard; @grimm2010odd]. The model is written in the Rust
[@rust] language, which easily provides fast and parallel execution, high-level
constructs, and type safety that prevents many types of programming errors. An
initial attempt to implement the model in Python [@python] suffered from very
slow execution even after using Cython [@cython] for some of the more
time-critical functions. Using literate programming [@knuth1984literate] to the
extent practical in Rust, the same file that can be compiled using the Rust
compiler also directly generates this model description.

## Purpose

The dispersal model generates the spatial distribution of cultures (together
with the history of these cultures) arising as emergent property from low-level
demographic and migratory processes. The demographic and migratory processes are
driven by the availability of resources, which is higher where culture-mediated
cooperation occurs. As such, there is a feedback loop between resource
availability, cultural proximity, and cooperation, which drives the population
dynamics of hunter-gatherers in the model.

In the current, first stage, the purpose of the model is to investigate patterns
languages disperse and split, driven only by the necessary interactions between
humans. The model is designed with extension to more concrete research questions
in mind. It is structured to be easily applied to study the history of the
settlement of the Americas at a later time, but in its current iteration assumes
constant climate and as such cannot be expected to produce directly comparable
results.

*/
// TODO: Dear Rustacean, I know that my use of documentation comments is
// hazardous, because they land in the generated documentation in a different
// order, or attached to things that are only accidentally coming afterwards.
// The structure is warranted by direct inclusion of the source in the paper
// text.

// Load useful modules
use std::f64::consts::PI;

use std::thread;
use rand::prelude::*;
use rayon::prelude::*;
use std::collections::HashMap;

mod util;
use util::{greater_of_two, smaller_of_two};

mod debug;
mod ecology;
mod tests;

use submodels::parameters::Parameters;

/**

## Entities, state variables, and scales

To construct a demographic migration model with culture, for the purpose of the
project we set out here, we need the following ingedients.

 - A representation of culture, which is at the least able to undergo neutral
   evolution (showing heritability and random variation, not necessarily
   fitness) and which in addition allows horizontal transfer of cultural traits
   other than from a parent to a child population.
 - Agents that carry cultural traits and are located in geographical space,
   which has differing ecological features
 - A system that drives the demographics of the agents in time and space
 - A way for culture and population dynamics to interact in a way that can
   create create distinct cultural areas instead of a vast cultural cline or
   dialect continuum.

Here we define those agents, cultures, and the underlying
ecology. The processes that drive the simulation will be explained in the
subsequent steps.

The model consists of agents interacting on a hexagonal discrete global grid in
discrete time. One time step is supposed to model a season with a duration of
half a year.

*/
type HalfYears = u32;
/**

### Families

The main decision-making agents of the simulation are families
[@crema2014simulation, I think?]. Families can migrate between grid cells and
form links to other families in the context of cooperation to extract resources.

1. Survival (and generation of off-splits) of a family and its culture depends
   on the available resources
2. Migration is resource-driven, and families migrate in their entirety
3. Families interact with each other to improve their chances of survival, which
   leads to an assimilation of cultural traits

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
    location: hexgrid::Index,
    /// The previous locations of the agent. This is useful for some bits of
    /// analysis, but agents are also guaranteed to have knowledge about their
    /// recent locations within range, where other knowledge depends on random
    /// factors.
    location_history: Vec<hexgrid::Index>,
    /// The amount of stored resources, in kcal, the family has access to.
    stored_resources: KCal,

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
their chances at extracting resources. (Following XXX, cooperation could also
mean sharing resources between different families. To simplify the model, this
implementation does not contain that effect.) The cooperative groups formed by
cooperating families are higher-level agents created ad-hoc in each time step.
They do not persist or have effect beyond a single time step. Resource
exploitation happens at the level of the cooperative group and is distributed
to the individual families after the fact.


### Cultures

Every family has a culture. These are very abstract and vastly simplified, due
to the lack of quantitative data on cultural evolution in a framework comparable
to the one used for this study. Based on the need to have a data structure that
supports random drift, a low but positive chance of back-mutation, and a large
number of equally similar cultures, we describe culture using a binary vector. A
detailed discussion of this choices and the alternatives can be found in
[#submodel:culture]. Computers make natural use of the binary representation to
store integer numbers, so a binary vector is equivalent to an unsigned integer
of sufficient size, which is faster to use in computations and more efficient to
store.
 
*/
type Culture = u64;
/**

### Grid and Patches

The geography of the simulation is described by a hexagonal equal-area discrete
global grid. The grid logic is implemented in Uber's H3 library [@uber-h3] and
exposed to the simulation by the `hexgrid` module.

*/
mod hexgrid;
/**

Each individual grid cell contains exactly one patch. Each of these patches has
a maximum and current availability of resources, measured in kcal available over
the course of half a year. In general, whereever possible, resources are
measured in kcal (in SI units: 1 kcal = 4.184 kJ)

*/
type KCal = f32;
// TODO: Maybe KCal should be a finite f32? If resources become infinite,
// something has gone wrong. I bet there is a Rust crate for that…
/**

In the current first iteration of the model, the resource maximum is constant
throughout the simulation. This is a natural point of extension in the future:
Future work might add seasonality (with every other time step corresponding to
the more plentiful half of the year and the others to the scarcer half) or
include paleoclimate data, but that only becomes relevant once the model shows
fundamentally reasonable dynamics.
 
*/
pub struct Patch {
    /// The resources currently available (not necessarily accessible, see XXX)
    /// in this patch, in kcal / season.
    resources: KCal,
    /// The maximum possible resources available in a time step, also in kcal /
    /// season.
    max_resources: KCal,
}
/**

### State

The associations between gridcells and patches and between gridcells and the
families located there (stored in the Family's `location`) are the core of the
model state. The state also tracks the time, measured in time steps
corresponding to half a year each, since the start of the simulation.
 
*/
struct State {
    /// The patches of the model, indexed by their address according to the H3
    /// discrete global grid system. This set is fixed, although the properties
    /// of the patches may change with time.
    patches: HashMap<hexgrid::Index, Patch>,
    /// The agents (families) currently active in the simulation. This vector
    /// changes over time as new families are added and old families die out.
    families: Vec<Family>,
    /// The current time step, in half years since start of the simulation.
    t: HalfYears,
}
// TODO: At some point, we thought about storing the model parameters in the state, such that serializing and deserializing the state is exactly what is necessary to store/resume the simulation. Do it, or why not?
/**

## Process overview and scheduling

> Who (i.e., what entity) does what, and in what order? When are state variables
> updated? How is time modeled, as discrete steps or as a continuum over which
> both continuous processes and discrete events can occur? Except for very
> simple schedules, one should use pseudo-code to describe the schedule in every
> detail, so that the model can be re-implemented from this code. Ideally, the
> pseudo-code corresponds fully to the actual code used in the program
> implementing the ABM.

The model progresses in discrete time steps, each corresponding to half a year
of simulated time. The entire simulation consists of repeating the step to
simulate 15000 years.

The structure of a single time step consist of two parts, as follows.

 
*/
fn step(
    families: &mut Vec<Family>,
    patches: &mut HashMap<hexgrid::Index, Patch>,
    p: &Parameters,
    t: HalfYears,
) -> Option<HashMap<hexgrid::Index, Vec<Family>>> {
    let mut families_by_location = step_part_1(families, patches, p, t)?;
    step_part_2(&mut families_by_location, patches, p);
    Some(families_by_location)
}
/**

The first part focuses on the individual families, which shrink, grow, die,
split, and move. It constructs the mapping of families at the end of the season,
grouped by their location after potential moves. Because the movement of a
family depends on the distribution of othe families at he start of the season,
it can happen entirely in parallel.

TODO: This is supposed to be readable [pseudo]code. Make it really really readable, and bigger.

*/
fn step_part_1(
    families: &mut Vec<Family>,
    patches: &HashMap<hexgrid::Index, Patch>,
    p: &Parameters,
    t: HalfYears,
) -> Option<HashMap<hexgrid::Index, Vec<Family>>> {
    let mut families_by_location: HashMap<hexgrid::Index, Vec<Family>> = HashMap::new();
    let mut cultures_by_location: HashMap<hexgrid::Index, HashMap<Culture, usize>> = HashMap::new();

    for family in families.iter_mut() {
        if submodels::family_lifecycle::use_resources_and_maybe_shrink(
            &mut family.effective_size,
            &mut family.stored_resources,
            p,
        ) {
            family.seasons_till_next_child = greater_of_two(family.seasons_till_next_child, 2)
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
        thread::spawn(|| {observation::print_gd_cd(l_c)});
        let l_c = cultures_by_location.clone();
        thread::spawn(|| {observation::print_population_by_location(l_c)});
    }

    let mut rng = rand::thread_rng();
    families.shuffle(&mut rng);

    for mut family in families.drain(..) {
        if submodels::family_lifecycle::is_moribund(&family) {
            family.effective_size = 0;
            // Stop tracking this family, they are dead.
            continue;
        }

        let nearby = hexgrid::nearby_locations(family.location);

        match submodels::family_lifecycle::maybe_procreate(&mut family) {
            None => {}
            // In terms of scheduling a new family can (and if possible
            // should) move immediately when created. This behaviour is
            // taken from del Castillo (2013). Also, a descendant family
            // will move directly before their progenitor.
            Some(descendant) => {
                let observed = sensing::observe_neighbors(
                    &descendant,
                    patches,
                    &cultures_by_location,
                    &nearby,
                    p,
                );
                let d = adaptation::decide_on_moving(&descendant, observed, true, p);
                let destination = d.unwrap_or(family.location);

                // Update cultures_by_location
                let old_count = cultures_by_location
                    .get_mut(&family.location)?
                    .get_mut(&family.culture)?;
                *old_count -= descendant.effective_size;
                let new_count = cultures_by_location
                    .entry(family.location)
                    .or_insert_with(HashMap::new)
                    .entry(family.culture)
                    .or_insert(0);
                *new_count += descendant.effective_size;

                families_by_location
                    .entry(destination)
                    .or_insert(vec![])
                    .push(descendant);
            }
        }
        let observed = sensing::observe_neighbors(
            &family,
            patches,
            &cultures_by_location,
            &nearby,
            p,
        );
        let d = adaptation::decide_on_moving(&family, observed, false, p);

        // Update cultures_by_location
        let old_count = cultures_by_location
            .get_mut(&family.location)?
            .get_mut(&family.culture)?;
        *old_count -= family.effective_size;
        let new_count = cultures_by_location
            .entry(family.location)
            .or_insert_with(HashMap::new)
            .entry(family.culture)
            .or_insert(0);
        *new_count += family.effective_size;

        let destination = d.unwrap_or(family.location);
        families_by_location
            .entry(destination)
            .or_insert(vec![])
            .push(family);
    }
    Some(families_by_location)
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
    families_by_location: &mut HashMap<hexgrid::Index, Vec<Family>>,
    patches: &mut HashMap<hexgrid::Index, Patch>,
    p: &Parameters,
) -> Option<()>{
    for (patch_id, families) in families_by_location {
        assert!(families.len() > 0);
        let mut patch = patches.entry(*patch_id).or_insert(Patch {
            // TODO Should we actually panic in this case? After all, some
            // family moved to a patch that doesn't exist, how did they even
            // find it?
            resources: 0.0,
            max_resources: 0.0,
        });
        let mut resource_reduction: KCal = 0.0;
        let (mut cooperatives, sum_labor) = collectives::cooperatives(families.iter_mut().collect(), p)?;

        for mut cooperating_families in cooperatives.iter_mut() {
            resource_reduction += submodels::ecology::extract_resources(
                &mut patch,
                &mut cooperating_families,
                sum_labor,
                p,
            );
        }

        for cooperating_families in cooperatives {
            let families = cooperating_families.families;
            crate::submodels::culture::adjust(families, p);
        }
        submodels::ecology::exploit(&mut patch, resource_reduction);
    }

    for mut patch in patches.values_mut() {
        submodels::ecology::recover(&mut patch, p.resource_recovery, p.accessible_resources);
    }
    Some(())
}

/**
## Design concepts

### Basic priciples

> Which general concepts, theories, hypotheses, or modeling approaches are
> underlying the model’s design?

According to its purpose of providing a model for language dispersal and split,
our model draws on existing publications looking at for cultures changing in
time and space, focussing on taking a bottom-up approach on languages splitting.
A major useful reference is the PSMED model [@delcastillo2013modeling; @barcelo2013psmed; @barcelo2015simulating], which includes two drivers of cultural
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
mod crema2014simulaton {}
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
### Emergence

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
        cultures_by_location: &HashMap<hexgrid::Index, HashMap<Culture, usize>>,
        max_geo_distance: i32,
    ) -> HashMap<i32, HashMap<u32, usize>> {
        let mut gd_cd: HashMap<i32, HashMap<u32, usize>> = HashMap::new();
        for (l1, cs1) in cultures_by_location {
            for (l2, cs2) in cultures_by_location {
                let gd = hexgrid::hex_distance(*l1, *l2);
                if gd > max_geo_distance {
                    continue;
                }
                for (c1, count1) in cs1 {
                    for (c2, count2) in cs2 {
                        if l1 == l2 && c1 == c2 {
                            let count = gd_cd
                                .entry(gd)
                                .or_insert_with(HashMap::new)
                                .entry(0)
                                .or_insert(0);
                            *count += count1 * (count2 + 1) / 2;
                            break;
                        }
                        let cd = submodels::culture::distance(*c1, *c2);
                        let count = gd_cd
                            .entry(gd)
                            .or_insert_with(HashMap::new)
                            .entry(cd)
                            .or_insert(0);
                        *count += count1 * count2;
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
    fn crema_fission_fusion_analysis() {
        // TODO: write this
    }

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
    
    pub fn similar_culture(c1: Culture, c2: Culture, p: &Parameters) -> bool {
        submodels::culture::distance(c1, c2) < p.cooperation_threshold
    }
}

/**
### Adaptation

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
mod adaptation {
    use crate::*;

    pub fn decide_on_moving<'a>(
        family: &'a Family,
        known_destinations: Vec<(hexgrid::Index, &Patch, usize, usize)>,
        avoid_stay: bool,
        p: &Parameters,
    ) -> Option<hexgrid::Index> {
        let mut kd = known_destinations.iter();
        let (mut target, mut patch, mut cooperators, mut competitors): (
            hexgrid::Index,
            &Patch,
            usize,
            usize,
        );
        {
            let (t, p, c, d) = kd.next()?;
            target = *t;
            patch = p;
            // This is the patch that contains the family, correct fo that.
            cooperators = *c - family.effective_size;
            competitors = *d;
        }
        if avoid_stay {
            let (t, p, c, d) = kd.next()?;
            target = *t;
            patch = p;
            cooperators = *c;
            competitors = *d;
        }
        let mut max_gain: KCal = prediction::expected_resources_from_patch(
            family.effective_size as f32,
            patch,
            cooperators as f32,
            competitors as f32,
            p,
        );

        let threshold = if avoid_stay {
            0.
        } else {
            max_gain + p.time_step_energy_use * p.evidence_needed
        };

        objectives::best_location(family.effective_size, &mut max_gain, kd, p, threshold)
            .or(Some(target))
    }
}

/**
### Objectives

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
expected resource gain. If multiple locations are equally good (up to
floating point accuracy), take one of those at random.

*/
    // That is, this function omputes the argmax of
    // `expected_resources_from_patch`, drawing at random between equal options.
    pub fn best_location(
        size: usize,
        max_gain: &mut KCal,
        kd: std::slice::Iter<(hexgrid::Index, &Patch, usize, usize)>,
        p: &Parameters,
        threshold: KCal,
    ) -> Option<hexgrid::Index> {
        let mut rng = rand::thread_rng();

        // This variable `c` is used to randomly draw between several
        // equally-optimal options. It counts the number of best options
        // encountered so far.
        let mut c = 0;
        let mut target: Option<hexgrid::Index> = None;
        for (coords, patch, cooperators, competitors) in kd {
            let expected_gain = prediction::expected_resources_from_patch(
                size as f32,
                patch,
                *cooperators as f32,
                *competitors as f32,
                p,
            );
            if expected_gain >= *max_gain {
                if expected_gain < threshold {
                    continue;
                }
                if (expected_gain - *max_gain).abs() < std::f32::EPSILON {
                    c += 1;
                    if rng.gen_range(0, c + 1) < c {
                        continue;
                    }
                } else {
                    c = 0
                }
                target = Some(*coords);
                *max_gain = expected_gain;
            }
        }
        target
    }
}
/**
### Learning

> Many individuals or agents (but also organizations and institutions) change
> their adaptive traits over time as a consequence of their experience? If so,
> how?

There is a minor effect of experience on location choice in that a family knows
nearby locations only with a small random chance, but they always have knowledge
about the locations they were at in the last 4 years (8 time steps).

TODO: Put that 4 in a static variable and put a doctest here that confirms it.
 
*/
mod learning {
    use crate::*;

    pub fn known_location(
        history: &Vec<hexgrid::Index>,
        nearby: &Vec<hexgrid::Index>,
        attention_probability: f32,
    ) -> Vec<hexgrid::Index> {
        let mut result = vec![];
        for location in nearby {
            if history[0..smaller_of_two(8, history.len())].contains(&location) {
                result.push(*location);
                continue;
            }
            if sensing::attention(attention_probability) {
                result.push(*location);
                continue;
            }
        }
        result
    }
}

/**
### Prediction

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
    use crate::*;

    pub fn expected_resources_from_patch(
        my_effective_size: f32,
        patch: &Patch,
        cooperators: f32,
        competitors: f32,
        p: &Parameters,
    ) -> KCal {
        assert!(my_effective_size + cooperators > 0.);
        submodels::ecology::resources_from_patch(
            patch,
            my_effective_size + cooperators,
            competitors,
            true,
            p,
        ) * my_effective_size / (my_effective_size + cooperators)
    }
}

/**
### Sensing

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
mod sensing {
    use crate::*;

    pub fn attention(attention_probability: f32) -> bool {
        random::<f32>() < attention_probability
    }

    pub fn scout<'a>(
        location: hexgrid::Index,
        reference_culture: Culture,
        patches: &'a HashMap<hexgrid::Index, Patch>,
        cultures_by_location: &HashMap<hexgrid::Index, HashMap<Culture, usize>>,
        p: &Parameters
    ) -> Option<(hexgrid::Index, &'a Patch, usize, usize)> {
        let cultures = match cultures_by_location.get(&location) {
            None => return Some((location, patches.get(&location)?, 0, 0)),
            Some(c) => c,
        };
        let mut cooper: usize = 0;
        let mut compet: usize = 0;
        for (culture, count) in cultures {
            if emergence::similar_culture(reference_culture, *culture, p) {
                cooper += count;
            } else {
                compet += count;
            }
        }
        Some((location, patches.get(&location)?, cooper, compet))
    }

    pub fn observe_neighbors<'a>(
        family: &'a Family,
        patches: &'a HashMap<hexgrid::Index, Patch>,
        cultures_by_location: &HashMap<hexgrid::Index, HashMap<Culture, usize>>,
        neighbors: &Vec<hexgrid::Index>,
        p: &Parameters
    ) -> Vec<(hexgrid::Index, &'a Patch, usize, usize)> {
        let mut result: Vec<(hexgrid::Index, &'a Patch, usize, usize)> = vec![];
        match scout(
            family.location,
            family.culture,
            patches,
            cultures_by_location,
            p
        ) {
            None => {}
            Some((location, patch, cooperators, competitors)) => {
                result.push((location, patch, cooperators, competitors))
            }
        }
        for location in
            learning::known_location(&family.location_history, neighbors, p.attention_probability)
        {
            if location == family.location {
                continue;
            }
            match scout(location, family.culture, patches, cultures_by_location, p) {
                None => {}
                Some((location, patch, cooperators, competitors)) => {
                    result.push((location, patch, cooperators, competitors));
                }
            }
        }
        result
    }
}

/**
### Interaction

> What kinds of interactions among agents are assumed? Are there direct
> interactions in which individuals encounter and affect others, or are
> interactions indirect, e.g., via competition for a medi- ating resource? If
> the interactions involve communication, how are such communications
> represented?

Agents compete for limited resources. There is a maximum of resources that can
be extracted from a patch, and the extracted resources are distributed among all
agents in that patch. This distribution is not completely even: Members of a
bigger group of cooperators competing with a smaller cooperative gets an
over-proportional part of the total extracted, measured by their effective size.
 
*/
mod interaction {
    use crate::*;
    pub fn split_resources(
        my_relative_returns: KCal,
        others_relative_returns: KCal,
        total_resources_available: KCal,
    ) -> KCal {
        assert!(my_relative_returns + others_relative_returns > 0.);
        let my_part = my_relative_returns / (my_relative_returns + others_relative_returns)
            * smaller_of_two(
                my_relative_returns + others_relative_returns,
                total_resources_available,
            );
        assert!(my_part.is_finite());
        my_part
    }
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
        1. + (n_cooperators - 1.).powf(cooperation_gain) / 10.
    }
}

/**
### Stochasticity

> What processes are modeled by assuming they are random or partly random? Is
> stochasticity used, for example, to reproduce variability in processes for
> which it is unimportant to model the actual causes of the variability? Is it
> used to cause model events or behaviors to occur with a specified frequency?

TODO: Is there a way to get rust to aggregate all functions that use it? Can I
maybe somehow abuse the trait system to all those locations show up here, with
comments?
 
*/
/**
### Collectives

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
        pub total_efficiency: f32,
    }

    pub fn cooperatives<'a>(
        families_in_this_location: Vec<&'a mut Family>,
        p: &Parameters
    ) -> Option<(Vec<Cooperative<'a>>, f32)> {
        let mut groups: Vec<Cooperative> = vec![];

        for family in families_in_this_location {
            let mut joined_group: Option<&mut Cooperative> = None;
            for group in groups.iter_mut() {
                let mut join = true;
                for other_family in group.families.iter() {
                    if !emergence::similar_culture(family.culture, other_family.culture, p) {
                        join = false;
                        break;
                    }
                }
                if join {
                    joined_group = Some(group);
                    break;
                }
            }
            let size = family.effective_size as f32;
            match joined_group {
                None => {
                    groups.push(Cooperative{
                        families: vec![family],
                        total_efficiency: size
                    });
                }
                Some(group) => {
                    group.families.push(family);
                    group.total_efficiency += size;
                }
            }
        }

        let mut rng = rand::thread_rng();
        use rand_distr::Normal;
        let dist = Normal::new(0., 1.).ok()?;
        let mut sum_labor = 0.;
        for group in groups.iter_mut() {
            group.total_efficiency = group.total_efficiency + dist.sample(&mut rng) * p.payoff_std / group.total_efficiency.powf(0.5);
            group.total_efficiency = greater_of_two(group.total_efficiency, 0.);
            sum_labor += group.total_efficiency;
        }
        Some((groups, sum_labor))
    }
}

/*
## Observation

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
        cultures_by_location: HashMap<hexgrid::Index, HashMap<Culture, usize>>,
    ) {
        println!(
            "POPULATION: {:?}",
            cultures_by_location
                .par_iter()
                .map(|(k, v)| {
                    let g = hexgrid::geo_coordinates(*k);
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
    pub fn print_gd_cd(cultures_by_location: HashMap<hexgrid::Index, HashMap<Culture, usize>>) {
        let gd_cd: HashMap<i32, HashMap<u32, usize>> =
            emergence::cultural_distance_by_geographical_distance(&cultures_by_location, 5 * 5 * 2);
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
fn initialization(precipitation: &Vec<u16>, width: usize, p: &Parameters) -> Option<State> {
    let start1: hexgrid::Index = hexgrid::closest_grid_point(-159.873, 65.613)?;
    let start2: hexgrid::Index = hexgrid::closest_grid_point(-158.2718, 60.8071)?;

    let area: f32 = libh3::hex_area_km_2(5) as f32;

    let mut patches: HashMap<hexgrid::Index, Option<Patch>> = HashMap::new();
    let mut new_patches: std::collections::BinaryHeap<hexgrid::Index> =
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
        let geo = hexgrid::geo_coordinates(next);
        let latitude = geo.lat * 180. / PI;
        let longitude = geo.lon * 180. / PI;
        if p.boundary_west < longitude
            && longitude < p.boundary_east
            && p.boundary_south < latitude
            && latitude < p.boundary_north
        {
            let patch = ecology::patch_from_coordinates(geo, precipitation, width);
            match patch {
                None => {
                    patches.insert(next, None);
                    continue;
                }
                Some(resources) => {
                    patches.insert(
                        next,
                        Some(Patch {
                            resources: resources
                                * p.time_step_energy_use
                                * area
                                * p.accessible_resources,
                            max_resources: resources * p.time_step_energy_use * area,
                        }),
                    );
                }
            };
            for q in hexgrid::nearby_locations(next) {
                if patches.contains_key(&q) {
                    continue;
                }
                new_patches.push(q);
            }
        } else {
            patches.insert(next, None);
        }
    }
    let f1 = Family {
        descendence: String::from("A"),
        location: start1,
        location_history: vec![],
        seasons_till_next_child: 4,
        culture: 0b000_000_000_000_000,

        effective_size: 2,
        number_offspring: 0,
        seasons_till_next_mutation: None,
        stored_resources: 1_600_000.,
    };
    let f2 = Family {
        descendence: String::from("F"),
        location: start2,
        location_history: vec![],
        seasons_till_next_child: 4,
        culture: 0b111_111_111_111_111,

        effective_size: 2,
        number_offspring: 0,
        seasons_till_next_mutation: None,
        stored_resources: 1_600_000.,
    };
    // let families: HashMap<hexgrid::Index, Vec<Family>> = HashMap::new();
    // families.insert(start1, vec![f1]);
    // families.insert(start2, vec![f2]);
    Some(State {
        patches: patches
            .drain()
            .filter_map(|(i, p)| match p {
                None => None,
                Some(q) => Some((i, q)),
            })
            .collect(),
        families: vec![f1, f2],
        t: 0,
    })
}

/**
## Input Data

> Does the model use input from external sources such as data files or other
> models to represent processes that change over time?
 
*/
mod input {
    // None at the moment. It might come, though, when we use paleoclimate data.
}

/**
## Submodels

> What, in detail, are the submodels that represent the processes listed in
> ‘Process overview and scheduling’? What are the model parameters, their
> dimensions, and reference values? How were submodels designed or chosen, and
> how were they parameterized and then tested?

*/
mod submodels {
    pub mod culture {
        use crate::{Culture, Family, Parameters};
        use rand::prelude::*;

        pub fn distance(c1: Culture, c2: Culture) -> u32 {
            (c1 ^ c2).count_ones()
        }

        pub fn mutate_culture(family: &mut Family, p: &Parameters) {
            match family.seasons_till_next_mutation {
                None => {
                    family.seasons_till_next_mutation =
                        Some(random::<f64>().log(1. - p.culture_mutation_rate) as u32);
                    mutate_culture(family, p);
                },
                Some(0) => {
                    let i: u8 = rand::thread_rng().gen_range(0, p.culture_dimensionality);
                    family.culture ^= 1 << i;
                    family.seasons_till_next_mutation =
                        Some(random::<f64>().log(1. - p.culture_mutation_rate) as u32);
                },
                Some(k) => {
                    family.seasons_till_next_mutation = Some(k-1);
                }
            }
        }

        pub fn adjust(mut cooperating_families: Vec<&mut Family>, p: &Parameters) {
            for f in &mut cooperating_families {
                mutate_culture(f, p);
            }
            let family: Option<&&mut Family> = cooperating_families.choose(&mut rand::thread_rng());
            match family {
                None => {}
                Some(f) => {
                    let target: Culture = f.culture;
                    for f2 in cooperating_families {
                        f2.culture = target;
                    }
                }
            }
        }
    }

    pub mod family_lifecycle {
        use crate::{Family, KCal, Parameters};

        pub fn resources_at_season_end(resources: KCal, size: usize, p: &Parameters) -> KCal {
            let mut resources_after: KCal = resources - (size as f32) * p.time_step_energy_use;
            if resources_after > 0. {
                resources_after *= 1. - p.storage_loss;
            }
            resources_after
        }

        pub fn maybe_procreate(family: &mut Family, ) -> Option<Family> {
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
        ) -> bool {
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
        use crate::{interaction, KCal, Parameters, Patch, collectives};

        pub fn resources_from_patch(
            patch: &Patch,
            labor: f32,
            others_labor: f32,
            _estimate: bool,
            p: &Parameters,
        ) -> KCal {
            let my_relative_returns: f32 = p.time_step_energy_use * labor * interaction::effective_labor_through_cooperation(labor, p.cooperation_gain);
            let others_relative_returns = if others_labor > 0. {
                p.time_step_energy_use * others_labor * interaction::effective_labor_through_cooperation(others_labor, p.cooperation_gain)
            } else {
                0.
            };
            interaction::split_resources(
                my_relative_returns,
                others_relative_returns,
                patch.resources,
            )
        }

        pub fn exploit(patch: &mut Patch, resource_reduction: KCal) {
            patch.resources -= resource_reduction;
            if patch.resources < 0. {
                patch.resources = 0.;
            }
        }

        pub fn recover(patch: &mut Patch, resource_recovery: f32, accessible_resources: f32) {
            if patch.resources < patch.max_resources * accessible_resources {
                let underlying_resources =
                    patch.resources + (1. - accessible_resources) * patch.max_resources;
                patch.resources += underlying_resources
                    * resource_recovery
                    * (1. - underlying_resources / patch.max_resources);
            }
            assert!(patch.resources.is_normal());
        }

        pub fn extract_resources(
            patch: &mut Patch,
            group: &mut collectives::Cooperative,
            total_labor_here: f32,
            p: &Parameters,
        ) -> KCal {
            let labor: f32 = group.total_efficiency;

            assert!(labor > 0.);
            let resources_extracted =
                resources_from_patch(&patch, labor, total_labor_here - labor, false, &p);
            for family in group.families.iter_mut() {
                family.stored_resources +=
                    resources_extracted * family.effective_size as f32 / labor as f32;
            }
            // This function really does not belong here
            resources_extracted
        }
    }

    pub mod parameters {
        use crate::KCal;

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

            pub boundary_west: f64,
            pub boundary_east: f64,
            pub boundary_south: f64,
            pub boundary_north: f64,
        }
    }
}


fn main() {
    let mut p = Parameters {
        attention_probability: 0.1,
        time_step_energy_use: 2263. as KCal * 365.242_2 / 2.,
        storage_loss: 0.33,
        resource_recovery: 0.20,
        culture_mutation_rate: 6e-3,
        culture_dimensionality: 20,
        cooperation_threshold: 6,
        cooperation_gain: 0.5,
        accessible_resources: 0.2,
        evidence_needed: 0.3,
        payoff_std: 0.1,

        boundary_west: -168.571_541,
        boundary_east: -34.535_395,
        boundary_south: -56.028_198,
        boundary_north: 74.52671,
    };
    let mut max_t: HalfYears = 20000;
    parse_args(&mut p, &mut max_t);
    run(p, max_t);
}

fn parse_args(p: &mut Parameters, max_t: &mut HalfYears) {
    let mut parser = argparse::ArgumentParser::new();
    parser.set_description("Run a dispersal simulation");
    // TODO: Attach an ArgumentParser to the parameters object, so parameters can be set from the CLI.
    parser.refer(&mut p.attention_probability)
        .add_option(&["--attention-probability"], argparse::Store, "attention probability");
    parser.refer(&mut p.culture_mutation_rate)
        .add_option(&["--culture-mutation-rate"], argparse::Store, "culture mutation rate");
    parser.refer(&mut p.culture_dimensionality)
        .add_option(&["--culture-dimensionality"], argparse::Store, "culture dimensionality");
    parser.refer(&mut p.cooperation_threshold)
        .add_option(&["--cooperation-threshold"], argparse::Store, "threshold under which cooperation happens");
    parser.refer(&mut p.cooperation_gain)
        .add_option(&["--cooperation-gain"], argparse::Store, "exponent in benefits from cooperation");
    parser.refer(max_t)
        .add_option(&["--steps"], argparse::Store, "number of half years to simulate");
    parser.parse_args_or_exit();
}

fn run(p: Parameters, max_t: HalfYears) -> Option<()> {
    println!("# Initialization ...");
    // FIXME: There is something messed up in lat/long -> image pixels. It works
    // currently, but the names point out that I misunderstood something, eg.
    // with the layout of the tiff in memory or similar.

    let (vec, width) = ecology::load_precipitation_tif()?;
    println!("Loaded precipitation data");
    let mut s: State = initialization(&vec, width as usize, &p)?;
    println!("Initialized");

    loop {
        let families_by_location = step(&mut s.families, &mut s.patches, &p, s.t)?;
        for (location, families) in families_by_location {
            for mut family in families {
                family.location_history.push(family.location);
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
            return None;
        }
    }
    Some(())
}
