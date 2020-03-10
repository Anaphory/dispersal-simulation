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

import gavin2017processbased
import binford2001constructing
from util import *

# 2. Entities, state variables, and scales
# ========================================
@attr.s
class Patch:
    """A patch of land with resources

    Resources can vary through time, up to a fixed maximum
    """
    resources: kcal = attr.ib() # Resources available available over a period of 6 months
    max_resources: kcal = attr.ib() # Maximum esources available available over a period of 6 months

@attr.s
class Cells():
    """A grid of cells

    By the default implementation in this module, every cell is a georeferenced
    hex cells in the Americas, each containing a Patch, representing an area of
    ca. 450_000_000 m²

    """
    patches: Mapping[Index, Patch] = attr.ib(factory=dict)
    neighbors: Callable[[Index], Sequence[Index]]
    neighbors_within_distance: Callable[[Index, meters], Sequence[Index]]


Culture = Tuple[int, int, int, int, int, int, int, int, int, int, int]

@attr.s
class Family:
    """A family group agent

    Families are the decision-making agent in our model. Families can migrate
    between cells and form links to other Families in the context of
    cooperation to extract or share resources.

    """
    descendence: str = attr.ib()
    # The agent's history of decendence, also serving as unique ID
    location: Index = attr.ib()
    # A Cells index
    culture: Culture = attr.ib()
    number_offspring: int = attr.ib(default=0)
    # The number of descendant families this family has given rise to so far
    effective_size: int = attr.ib(default=2)
    # The effective size of the family in number of adults. One adult is
    # assumed to consume the same amount of food/energy and to contribute the same labor to
    # foraging as any other adult.
    stored_resources: kcal = attr.ib(default=0)
    # The amount of stored resources, in kcal, the family has access to
    plenty_seasons_since_last_child: int = attr.ib(0)
    # The number of plentiful seasons since the birth of the previous child

@attr.s
class State:
    """A model state"""
    grid: Cells = attr.ib()
    # The state of cells at this time step
    families: List[Family] = attr.ib(factory=list)
    # A list of families under simulation
    t: halfyears = attr.ib(default=0)
    # Time steps since simulation start, in 1/2 of a year.

class CooperativeGroup(list):
    """A group of Families cooperating"""
    efficiency: float = attr.ib() # The group's efficiency at exploiting resources, between 0 and 1

# 3. Process overview and scheduling
# ==================================
def step(state: State) -> State:
    """Run a simulation step."""
    inhabited_patches = DefaultDict[Index, List[Family]](list)
    # Keeping track of inhabited patches locally is easier than looking for them post-hoc.

    for family in shuffle(state.families):
        use_resources_and_maybe_shrink(family)
        maybe_grow(family)
        if is_moribund(family):
            state.families.remove(family)
            continue
        descendant = maybe_procreate(family)
        if descendant is not None:
            state.families.append(descendant)
            # In terms of scheduling, this means in particular that a new
            # family can move immediately. This behaviour is taken from del
            # Castillo (2013). Also, a descendant family will move *before*
            # their progenitor.
            decide_on_moving_and_maybe_move(
                descendant,
                state.grid.patches[family.location],
                observe_neighbors(
                    family, state.grid.patches, state.families,
                    state.grid.neighbors_within_distance))
            inhabited_patches[descendant.location].append(descendant)
        # Subsequent earlier moves affect the possible targets of later moves,
        # so scheduling matters and shuffling is important to remove
        # first-mover effects.

        decide_on_moving_and_maybe_move(
            family, state.grid.patches[family.location],
            observe_neighbors(
                family, state.grid.patches, state.families,
                state.grid.neighbors_within_distance))
        inhabited_patches[family.location].append(family)

    for patch_id, families in inhabited_patches.items():
        patch = state.grid.patches[patch_id]
        resource_reduction = 0
        # For exploitation, all parts of it happen at the same time.
        groups_and_efficiencies, sum_labor = cooperate(families)
        for cooperating_families in groups_and_efficiencies:
            resource_reduction += extract_resources(patch, cooperating_families, sum_labor)
            adjust_culture(cooperating_families)

        exploit(patch, resource_reduction)

    for patch in state.grid.patches.values():
        recover(patch)

    observation(state)

    state.t += 1

    return state

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
 • The family is generally the basic unit of decision-making and migration in
   hunter-gatherer societis

4.2 Emergence
=============

"""

params: Namespace

def geographical_distance(index1: Index, index2: Index) -> meters:
    m1, i1, j1 = index1
    m2, i2, j2 = index2
    return numpy.asarray(GEODESIC.inverse(
        coordinates[m1][i1, j1], coordinates[m2][i2, j2])[:, 0])[0]

def plot_cultural_distance_by_geographical_distance(state: State) -> Tuple[Sequence, Sequence]:
    """Plot the cultural distance vs. geographical distance

    We expect that the interplay of migration, cooperation and cultural
    similarity leads to regions of similar cultures with steep geographical
    transition areas to regions containing agents with different culture
    vectors.

    This means that the plot of cultural distances vs. geographical distances
    should therefore show small cultural distances for small geographical
    distances. There should be a critical geographical distance of cohesion
    where the distribution of cultural distances becomes bimodal, with one mode
    being lower than the cooperation threshold and one mode above the
    cooperation threshold. For large geographical distances, the cultural
    distances should be high, but with a big variance.

    """
    cult_dists: List[float] = []
    geo_dists = []
    for family1, family2 in itertools.combinations(
            state.families, 2):
        geo_dist = geographical_distance(family1.location, family2.location)
        cult_dist = cultural_distance(family1.culture, family2.culture)
        geo_dists.append(geo_dist); cult_dists.append(cult_dist)
    return geo_dists, cult_dists


"""
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
   remaining agents in their old position before moving. Agents extrapolate
   their resource availabiltiy with and without expending labour from the
   current state of the system, ie. without taking potential family growth or
   shrinking into account and without taking the movement of agents into
   account.
"""

def family_expects_current_patch_will_be_enough(family, resource_gain):
    return resources_at_season_end(
        family.stored_resources + resource_gain,
        family.effective_size) > 0
"""
4.7 Sensing
===========
 • Agents have access to the current state of their neighborhood (actual
   available resources and position, size and culture vectors of nearby agents)
   for the decision process whether and where to move.
   FIXME: Currently, the function that describes the agent decision process
   does not make knowledge about other agents available to the decider, due to
   current limitations of the architecture…
"""

def observe_neighbors(
        family: Family,
        patches: Mapping[Index, Patch],
        all_families: Sequence[Family],
        neighbor_generator: Callable) -> Iterator[
        Tuple[Index, Patch, int, int]]:
    for dest in neighbor_generator(family.location, meters(148413)):
        cooperators, competitors = 0, 0
        for f in all_families:
            if f.location == dest:
                if cultural_distance(f.culture, family.culture) < params.cooperation_threshold:
                    cooperators += f.effective_size
                else:
                    competitors += f.effective_size
        yield dest, patches[dest], cooperators, competitors

"""
4.8 Interaction
===============
 • Agents compete for limited resources
 • Agents of similar culture cooperate in the exploitation of resources
"""

def cultural_distance(c1: Culture, c2: Culture) -> float:
    return sum(abs(e1-e2) for e1, e2 in zip(c1, c2))

"""
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

def observation(
        state: State,
        to:IO[str]=sys.stdout,
        extensive:IO[str]=open("log", "w")) -> None:
    # Number of families (agents)
    report = {}
    extreport = {}
    report["t"] = state.t
    report["Number of families"] = len(state.families)
    if len(state.families) == 0:
        raise StopIteration
    report["Total population count"] = sum([family.effective_size for family in state.families])
    report["Median stored resources"] = numpy.median([family.stored_resources for family in state.families])
    extreport["Locations"] = {family.descendence: family.location
                           for family in state.families}
    print(report, file=to)
    print(extreport, file=extensive)

# 5. Initialization
# =================
# In addition to describing the initialization, we list the parameters of the model
# here. Parameters are read from arguments specified on the command line, with
# the following default values. The sources for the default values are given.

coordinates = gavin2017processbased.hexagonal_earth_grid(
    gavin2017processbased.americas,
    gavin2017processbased.area)


def patch_from_grid_index(index: Index) -> Patch:
    m, i, j = index
    longitude, latitude = coordinates[m][i,j]
    data = get_data(longitude, latitude)
    resources = binford2001constructing.TERMD2(**data) * gavin2017processbased.area / 1_000_000 * params.time_step_energy_use
    return Patch(resources, resources)


def parse_args(args: Sequence[str]) -> Tuple[halfyears, kcal]:
    parser = argparse.ArgumentParser(description="Run dispersal model")
    parser.add_argument(
        "--n-steps", type=halfyears, default=30_000,
        help="The number of half-year steps to run the simulation for. The "
        "default is 30000 half-years, because the current scientific "
        "consensus is that the settlement of the Americas began some 15000 "
        "years BP.") # Cite: 
    parser.add_argument(
        "--daily-energy", type=kcal, default=2263,
        help="Energy consumption per adult per day. According to Pontzer et "
        "al. (2012), the mean daily total energy expenditure among Hadza was "
        "2263 kcal (mean of men and women). Smith & Smith (2003) work with a "
        "staple diet containing 2390 kcal/day.")
    parser.add_argument(
        "--payoff-standarddeviation", type=int, default=0.1,
        help="The standard deviation of the foraging contribution distribution, relative to the expected foraging contribution of an individual.")
    parser.add_argument(
        "--cooperation-gain", type=int, default=0.5,
        help="The exponent of the additional efficiency of a group working together. If, eg., cooperation_gain==0.3, then two individuals working together contribute as much as 4 working separately, 3 as 7.2, 4 as 11, 5 as 15 etc.")
    parser.add_argument(
        "--storage-loss", type=float, default=0.33,
        help="The proportion of stored resources lost per time step to "
        "various random effects, such as pests, rot, spoilage, or leaving "
        "them behind when migrating. Morgan (2012) pulls a number of 33% per "
        "year (focussing on the winter season when it matters) out of thin air.")
    parser.add_argument(
        "--cooperation_threshold", type=float, default=6,
        help="The minimum cultural distance where cooperation breaks down")
    parser.add_argument(
        "--resource-recovery", type=float, default=1.1,
        help="The growth rate of a path's resources over half a year")

    global params
    params = parser.parse_args(args)
    params.time_step_energy_use = params.daily_energy * 365.24219 / 2

    return params.n_steps, params.daily_energy

def initialization() -> State:
    grid = Cells()
    grid.patches = OnDemandDict(patch_from_grid_index)
    grid.neighbors_within_distance = cached_generator(neighbors_within(coordinates))

    return State(grid=grid, families=[
        Family(
            descendence="A",
            location=(0, 28, 101),# Start around Fairbanks
            culture=(2,2,2,2,2,2,2,2,2,2,2),
            stored_resources=16000000),
        Family(
            descendence="B",
            location=(0, 28, 101),# Start around Fairbanks
            culture=(2,2,2,2,2,2,2,2,2,2,2),
            stored_resources=16000000)])

# 6. Input Data
# =============

# This model does not currently use input data to represent time-varying
# processes. At a later stage, the inclusion of paleoclimate data for the
# Americas is intended as input data.

# 7. Submodels
# ============
#

def resources_at_season_end(resources: kcal, size: int) -> kcal:
    """7.1 Resource use model

    Use up a family's resources, and modify the remainder due to storage loss.

    Every adult consumes `daily_energy` kcal of food per day, for the period of
    one time step, i.e. half a year.

    """
    resources_after = resources - (
        size * params.time_step_energy_use)
    if resources_after > 0:
        resources_after = (1 - params.storage_loss)
    return resources_after


def extract_resources(patch: Patch, group: CooperativeGroup, total_labor_here: int):
    labor = sum([family.effective_size
                 for family in group])
    resources_extracted = resources_from_patch(
        patch, labor, total_labor_here - labor)
    for family in group:
        family.stored_resources += resources_extracted * family.effective_size / labor
    return resources_extracted


def use_resources_and_maybe_shrink(family):
    resources, size = family.stored_resources, family.effective_size
    while resources_at_season_end(resources, size) < 0 and size > 0:
        size -= 1
        family.plenty_seasons_since_last_child = 0
    family.effective_size = size
    family.stored_resources = resources_at_season_end(
        family.stored_resources, family.effective_size)


def maybe_grow(family: Family) -> None:
    if family.plenty_seasons_since_last_child > 1:
        family.effective_size += 1
        family.plenty_seasons_since_last_child = 0
    else:
        family.plenty_seasons_since_last_child += 1


def is_moribund(family: Family) -> bool:
    if family.effective_size < 2:
        return True
    else:
        return False


def maybe_procreate(family: Family) -> Optional[Family]:
    if family.effective_size < 10:
        return None
    else:
        family.effective_size -= 2
        family.number_offspring += 1
        return Family(
            descendence="{:s}:{:d}".format(family.descendence, family.number_offspring),
            culture=family.culture,
            location=family.location)


def decide_on_moving_and_maybe_move(
        family: Family,
        current_patch: Patch,
        neighbor_cells: Iterator[Tuple[Index, Patch, int, int]]) -> None:
    target, patch, cooperators, competitors = next(neighbor_cells)
    assert patch == current_patch
    max_gain: kcal = resources_from_patch(
        current_patch, cooperators, competitors) * family.effective_size / cooperators
    if family_expects_current_patch_will_be_enough(family, max_gain):
        move = "leisure" if numpy.random.random() < 0.0005 else False
    else:
        move = "desperate"

    if not move:
        return
    else:
        max_gain = 0 # FIXME why do I have to set this? It *should be* that a
                     # population that grows to big needs to much of local
                     # resources, so that it should be better to move.
        for coords, patch, cooperators, competitors in neighbor_cells:
            expected_gain = resources_from_patch(
                    patch, family.effective_size + cooperators, competitors) * (
                        family.effective_size / (family.effective_size + cooperators))
            if target and expected_gain > max_gain:
                target = coords
                max_gain = expected_gain
        family.location = target


def cooperate(families: Sequence[Family]) -> Tuple[Sequence[CooperativeGroup], int]:
    cooperative_groups: List[CooperativeGroup] = []
    sum_labor = 0
    for family in families:
        sum_labor += family.effective_size
        for group in cooperative_groups:
            for other_family in group:
                if cultural_distance(
                        family.culture, other_family.culture) > params.cooperation_threshold:
                    break
            else:
                group.append(family)
                break
        else:
            cooperative_groups.append(CooperativeGroup([family]))
    for group in cooperative_groups:
        group.efficiency = 1/len(families)
    return cooperative_groups, sum_labor


def effective_labour_after_cooperation(labor):
    """Effective total labor contribution of a group of cooperators"""
    # From crema2014simulation, adapted
    return labor * (1 + (labor - 1) ** params.cooperation_gain)


def resources_from_patch(
        patch: Patch,
        labor: int,
        others_labor: int) -> kcal:
    # From crema2014simulation
    my_relative_returns = numpy.random.normal(
        loc=params.time_step_energy_use * effective_labour_after_cooperation(labor),
        scale=params.payoff_standarddeviation * params.time_step_energy_use / labor ** 0.5)
    if others_labor:
        others_relative_returns = numpy.random.normal(
            loc=params.time_step_energy_use * effective_labour_after_cooperation(others_labor),
            scale=params.payoff_standarddeviation * params.time_step_energy_use / others_labor ** 0.5)
    else:
        others_relative_returns = 0
    return numpy.minimum(
        my_relative_returns + others_relative_returns,
        patch.resources) * (my_relative_returns) / (my_relative_returns + others_relative_returns)

def adjust_culture(cooperating_families: CooperativeGroup) -> None:
    ...

def exploit(patch: Patch, resource_reduction: kcal) -> None:
    patch.resources -= resource_reduction
    if patch.resources < -1e-2:
        # They should be bigger than 0, but there may be rounding errors
        raise AssertionError("Patch was over-exploited")
    elif patch.resources < 0:
        patch.resources = 0


def recover(patch: Patch) -> None:
    if patch.resources < patch.max_resources - 1:
        patch.resources += (
            patch.resources *
            params.resource_recovery *
            (1 - patch.resources / patch.max_resources))

# Run the simulation
if __name__ == "__main__":
    parse_args(["--n", "30000"])
s = initialization()
simulate(s, params.n_steps)

import cartopy.crs as ccrs
plt.gcf().set_size_inches(15, 15)
ax = plt.axes(projection=ccrs.PlateCarree())
ax.coastlines("50m")
ax.set_extent(gavin2017processbased.americas)


coords = []
for family in s.families:
    m, i, j = family.location
    coords.append(coordinates[m][i, j])
import matplotlib.pyplot as plt
plt.scatter(*zip(*coords), alpha=0.1)
plt.show()

sources = """
Bibliography
============

Grimm, Volker & Berger, Uta & Bastiansen, Finn & Eliassen, Sigrunn & Ginot,
    Vincent & Giske, Jarl & Goss-Custard, John et al. 2006. A standard protocol
    for describing individual-based and agent-based models. Ecological
    Modelling 198(1). 115–126. (doi:10.1016/j.ecolmodel.2006.04.023)

Grimm, Volker & Berger, Uta & DeAngelis, Donald L. & Polhill, J. Gary & Giske,
    Jarl & Railsback, Steven F. 2010. The ODD protocol: A review and first
    update. Ecological Modelling 221(23). 2760–2768.
    (doi:10.1016/j.ecolmodel.2010.08.019)

""".split("\n\n")[1:]
