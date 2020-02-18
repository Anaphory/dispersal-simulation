"""settlement

Model the settlement of the Americas using an agent-based model on a hex grid
with resource depletion, migration, and culture-dependent cooperation.

Model description following the ODD (Overview, Design concept, Details)
protocol (Grimm et al., 2006; Grimm et al., 2010) inline throughout the source
code of the model.

1. Purpose
==========

The settlement model generates a deep phylogeny of hunter-gatherer cultures
based on culture-mediated cooperation and resource-driven migration. The
summary statistics of this phylogeny (in particular diversification rates) are
to be compared to those known from language evolution. The model is structured
to be easily transferable to the history of the settlement of the Americas, but
requires paleoclimate data to produce results that can be compared to that
history.

"""
# General imports, which in Python should come on top of the module to make
# dependencies transparent.

import sys
import attr
import numpy
import pickle
from util import *

# Import some modules for the Entities, in particular types used to describe
# the Entities
Index = Tuple
meters = float
kcal = float
halfyears = int

# 2. Entities, state variables, and scales
# ========================================
@attr.s
class Patch:
    """A patch of land with resources

    """
    resources: kcal = attr.ib() # Resources available available over a period of 6 months
    max_resources: kcal = attr.ib() # Maximum esources available available over a period of 6 months

@attr.s
class Cells:
    """A grid of cells

    A topology of georeferenced hex cells in the Americas, each containing a Patch.

    By the default implementation in this module, every cell is a hex cell with an area of ca. 450_000_000 m²
    """
    patches: Mapping[Index, Patch] = attr.ib(factory=dict)
    neighbors: Callable[[Index], Sequence[Index]]
    neighbors_within_distance: Callable[[Index, meters], Sequence[Index]]

@attr.s
class Family:
    """A family group agent

    A family is decision-making agent in our model. Families can migrate
    between cells and form links to other Families in the context of
    cooperation to extract or share resources.

    """
    descendence: str = attr.ib()               # The agent's history of decendence, also serving as unique ID
    location: Index = attr.ib()                # A Cells index
    number_offspring: int = attr.ib(default=0) # The number of descendant families this family has given rise to so far
    effective_size: int = attr.ib(default=2)   # The effective size of the family in number of adults. One adult is assumed to consume 2000 kcal per day and to contribute the same labor to foraging as any other adult.

@attr.s
class State:
    """A model state"""
    grid: Cells = attr.ib()
    families: Sequence[Family] = attr.ib(factory=list)
    t: halfyears = attr.ib(default=0) # Time step, in 1/2 of a year.

@attr.s
class CooperativeGroup(list):
    """A group of Families cooperating"""
    efficiency: float = attr.ib() # The group's efficiency at exploiting resources, between 0 and 1

# 3. Process overview and scheduling
# ==================================
def step(state: State) -> State:
    """ Run a simulation step.

    """
    inhabited_patches = DefaultDict(list)
    for family in shuffle(state.families):
        spend_stored_resources(family)
        maybe_grow_or_shrink(family)
        if is_moribund(family):
            state.families.remove(family)
        # Subsequent earlier moves affect the possible targets of later moves,
        # so scheduling matters and shuffling is important to remove
        # first-mover effects.
        decide_on_moving_and_maybe_move(
            family,
            state.cells.neighbors_within_distance(family.location, 148413))
        inhabited_patches[family.location].append(family)

    for index, families in inhabited_patches:
        patch = state.grid.patches[index]
        resource_reduction = 0
        # For exploitation, all parts of it happen at the same time.
        groups_and_efficiencies = cooperate(families)
        for cooperating_families in groups_and_efficiencies:
            resource_reduction += extract_resources(
                patch, cooperating_families)
            adjust_culture(cooperating_families)

        exploit(patch, resource_reduction)

    for patch in state.grid.patches.values():
        recover(patch)

    observation(state)

    state.t += 1

def simulate(state: State, n_steps: int=15000 * 2) -> State:
    for t in range(n_steps):
        state = step(state)
    return state

"""
4. Design concepts
==================

This model draws heavily on two prior models.

4.1 Basic priciples
===================
 • Migration happens due to the need to find new resources in marginal evironments
 • Shared culture mediates cooperation, and cooperation fosters cultural cohesion
 • Using real(istic) geography provides better insight than abstract environments

4.2 Emergence
=============
 • We expect that the interplay of migration, cooperation and cultural
   similarity leads to regions of similar cultures with steep geographical
   transition areas to regions containing agents with different culture
   vectors.

4.3 Adaptation
==============
 • The adaptive behaviour of agents is restricted to their movement, which is
   implemented to optimize resource availability in the upcoming time step.

4.4 Objectives
==============
 • Agents evaluate the expected resources they will be able exploit in the
   current time step based on the resources available in patches within a
   certain distance from themselves, and the number of agents in each of those
   patches.
   FIXME: Currently, the function that describes the agent decision process
   does not make knowledge about other agents available to the decider, due to
   current limitations of the architecture…

4.5 Learning
============
 • Agents do not modify behaviour based on past experience.

4.6 Prediction
==============
 • Agents use the state of their neighborhood as proxy for the future state of
   the system for predicting expected number of resources available to them.
   Due to scheduling, this means that agents see on average half of the
   remaining agents in their new positions after moving and half of the
   remaining agents in their old position before moving.

4.7 Sensing
===========
 • Agents have access to the current state of their neighborhood (actual
   available resources and position, size and culture vectors of nearby agents)
   for the decision process whether and where to move.
   FIXME: Currently, the function that describes the agent decision process
   does not make knowledge about other agents available to the decider, due to
   current limitations of the architecture…

4.8 Interaction
===============
 • Agents compete for limited resources
 • Agents of similar culture cooperate in the exploitation of resources
 • Agents may avoid other agents when the expected gain from moving to their
   spot would be reduced due to their presence (different culture and thus
   competition instead of cooperation)

4.9 Collectives
===============
 • Groups of agents with similar culture on the same patch form temporary
   groups to collaborate for the exploitation of resources. The efficiency of
   the group is merely an emergent property of the individual Agent sizes.
   Adaptation of culture vectors is a function of these groups.

4.10 Observation
================
 • The following statistics are aggregated each time step.
"""

def observation(state: State, to: open=sys.stdout) -> None:
    # Number of families (agents)
    print(len(state.families), file=to)
    print(sum([family.effective_size for family in state.families]), file=to)

# 5. Initialization
# =================

def initialization() -> State:
    grid = Cells()
    return State(grid)

# 6. Input Data
# =============

# This model does not currently use input data to represent time-varying
# processes. At a later stage, the inclusion of paleoclimate data for the
# Americas is intended as input data.

# 7. Submodels
# ============

def spend_stored_resources(family: Family) -> None:
    ...

def maybe_grow_or_shrink(family: Family) -> None:
    ...

def is_moribund(family: Family) -> bool:
    return False

def decide_on_moving_and_maybe_move(family: Family, neighbor_cells: Sequence) -> None:
    ...

def cooperate(families: Sequence[Family]) -> Sequence[CooperativeGroup]:
    return CooperativeGroup(families)

def extract_resources(patch: Patch, cooperating_families: CooperativeGroup) -> kcal:
    return 0

def adjust_culture(cooperating_families: CooperativeGroup) -> None:
    ...

def exploit(patch: Patch, resource_reduction: kcal) -> None:
    ...

def recover(patch: Patch) -> None:
    ...


sources = """
Bibliography
============

Grimm, Volker, Uta Berger, Finn Bastiansen, Sigrunn Eliassen, Vincent Ginot,
    Jarl Giske, John Goss-Custard, et al. 2006. A standard protocol for
    describing individual-based and agent-based models. Ecological Modelling
    198(1). 115–126. doi:10.1016/j.ecolmodel.2006.04.023.

Grimm, Volker, Uta Berger, Donald L. DeAngelis, J. Gary Polhill, Jarl Giske &
    Steven F. Railsback. 2010. The ODD protocol: A review and first update.
    Ecological Modelling 221(23). 2760–2768.
    doi:10.1016/j.ecolmodel.2010.08.019.

""".split("\n\n")[1:]

