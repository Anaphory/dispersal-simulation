"""crema2014simulation

A replication of the fission-fusion model for hunter-gatherer population
dynamics and migration described by Crema (2014), with the option to add the
resource depletion described in scenario 4 described by Crema (2015).

Crema, Enrico R. 2014. A Simulation Model of Fission–Fusion Dynamics and
Long-Term Settlement Change. Journal of Archaeological Method and Theory 21(2).
385–404. (doi:10.1007/s10816-013-9185-4)

Crema, Enrico R. 2015. Modelling settlement rank-size fluctuations. In Wurzer,
Gabriel & Kowarik, Kerstin & Reschreiter, Hans (eds.), Agent-based Modeling and
Simulation in Archaeology (Advances in Geographic Information Science),
161–181. Cham: Springer International Publishing.
(doi:10.1007/978-3-319-00008-4_8) (https://doi.org/10.1007/978-3-319-00008-4_8)
(Accessed January 25, 2020.)

"""
import json
import numpy
import itertools
from matplotlib import pyplot as plt

model_parameters = dict(
    initial_agents = [10],
    time_steps = [500],
    patches = [100],
    basic_individual_payoff = [10],
    payoff_standard_deviation = [1.],
    # The paper calls this parameter ε “payoff variance”, but also describes
    # the process as “random draw from a Gaussian probability distribution with
    # […] a standard deviation ε”. Because the value is 1, it does not actually
    # make any difference.
    cooperative_benefit = [0.3, 0.5, 0.8],
    resource_input = [200],
    basic_reproduction_rate = [0.05],
    death_parameter_one = [0.8, 1.0, 1.2, 1.4],
    death_parameter_two = [5],
    spatial_interaction_range = [1, 100],
    decision_making_probability = [0.1, 0.5, 1.0],
    evidence_threshold = [3],
    observed_agents = [1e-7, 0.5, 1],
    resource_intrinsic_growth_rate = [2],
    resource_resilience = [None, 0.35, 0.4])


def fission(cell, neighbors, modifyable_population_grid):
    """Try to separate an individual from the group in the cell.

    Modify modifyable_population_grid (as a side effect), removing, if
    possible, one individual from `cell` and adding it to a random empty
    neighbor. If none of the neighbors are empty, do nothing.

    Examples
    ========
    >>> population = [0, 2, 0, 3, 2, 1, 0, 2]
    >>> fission(4, [3, 5], population)
    >>> population
    [0, 2, 0, 3, 2, 1, 0, 2]
    >>> fission(7, [6], population)
    >>> population
    [0, 2, 0, 3, 2, 1, 1, 1]
    >>> fission(1, [0, 2], population)
    >>> population[3:]
    [3, 2, 1, 1, 1]
    >>> population[1]
    1
    >>> population[0], population[2] in [(0, 1), (1, 0)]
    True

    Parameters
    ==========
    cell: (int, int)
        Or, more generally, an index into modifyable_population_grid.
        The index of the cell from which an indivual wants to fission.
    neighbors: list((int, int))
        Or, more generally, a sequence of indices.
        The neighboring cells, from which an empty one is drawn at random
        for the the individual to move to
    modifyable_population_grid: 2d numpy.array
       Or, more generally, a modifyable container
       The population counts of different cells

    Returns
    =======
    None

    Side effects
    ============
    modifies modifyable_population_grid, keeping the sum constant.
    """
    empty_cells = [
        n for n in neighbors
        if not modifyable_population_grid[n]]
    if not empty_cells:
        return
    modifyable_population_grid[cell] -= 1
    target = empty_cells[numpy.random.randint(len(empty_cells))]
    assert modifyable_population_grid[target] == 0
    modifyable_population_grid[target] = 1


def neighbors_within(cell, distance, shape, wrap=(False, False)):
    """Return all neighbors with a given Euclidean distance of cell.

    Examples
    ========
    >>> sorted(neighbors_within((1, 1), 1, (10, 10)))
    [(0, 1), (1, 0), (1, 2), (2, 1)]
    >>> sorted(neighbors_within((0, 0), 1, (10, 10)))
    [(0, 1), (1, 0)]
    >>> sorted(neighbors_within((0, 0), 1.5, (10, 10), (True, False)))
    [(0, 1), (1, 0), (1, 1), (9, 0), (9, 1)
    >>> sorted(neighbors_within((4, 4), 1.8, (10, 10), (True, False)))
    [(3, 3), (3, 4), (3, 5), (4, 3), (4, 4), (4, 5), (5, 3), (5, 4), (5, 5)]

    Parameters
    ==========
    cell: (int, int)
        A grid cell index
    distance: float
        The maximum Euclidean distance to still be included in the neighborhood
    shape: (int, int)
        The shape of the underlying array
    wrap: (bool, bool)
        whether to wrap the neighborhood around the borders. (False, False)
        gives a square, (True, True) a torus topology, everything else a
        non-twisted strip.
    """
    distance_sq = distance ** 2
    c1, c2 = cell
    for n1 in range(c1 - int(distance), c1 + int(distance) + 1):
        if not wrap[0] and (n1 < 0 or n1 >= shape[0]):
            continue
        for n2 in range(c2 - int(distance), c2 + int(distance) + 1):
            if not wrap[1] and (n2 < 0 or n2 >= shape[1]):
                continue
            if 0 < (n1 - c1) ** 2 + (n2 - c2) ** 2 <= distance_sq:
                yield n1 % shape[0], n2 % shape[1]


def simulation(
        initial_agents = 10,
        time_steps = 500,
        patches = 100,
        basic_individual_payoff = 10,
        payoff_standard_deviation = 1.,
        cooperative_benefit = 0.3,
        resource_input = 200,
        basic_reproduction_rate = 0.05,
        death_parameter_one = 0.8,
        death_parameter_two = 5,
        neighborhood = lambda x: neighbors_within(
            x, 1.5, (10, 10), (True, True)),
        spatial_interaction_range = 1.5,
        decision_making_probability = 0.1,
        evidence_threshold = 3,
        observed_agents = 1e-7,
        resource_resilience = 0.3,
        resource_intrinsic_growth_rate = 2,
        results = None):
    """Run Crema's simulation

    Returns
    =======
    list(numpy.array)
        The list of the state of the population in each time step
    """
    if resource_resilience is None:
        resource_threshold = 1.0
        change_resources = False
    else:
        resource_threshold = 1.0 - resource_resilience
        change_resources = True

    grid_shape = int(patches ** 0.5), patches // int(patches ** 0.5)
    population_grid = numpy.zeros(grid_shape, dtype=int)
    max_resource_grid = resource_input * numpy.ones_like(population_grid, dtype=float)
    resource_grid = max_resource_grid + 0
    for _ in range(initial_agents):
        population_grid[numpy.random.randint(population_grid.shape[0]),
                        numpy.random.randint(population_grid.shape[1])] += 1
    run_results = []
    for i in range(time_steps):
        foraging_contribution = numpy.random.normal(
            loc=population_grid * (
                basic_individual_payoff + (population_grid - 1)**cooperative_benefit),
            scale=payoff_standard_deviation / population_grid ** 0.5)
        fitness_grid = numpy.minimum(
            foraging_contribution,
            resource_grid * resource_threshold) / population_grid

        if change_resources:
            resource_grid -= fitness_grid
            resource_grid += (resource_intrinsic_growth_rate * resource_grid *
                              (1 - resource_grid / max_resource_grid))

        reproduction_prob = basic_reproduction_rate * (
            fitness_grid / basic_individual_payoff)
        numpy.nan_to_num(reproduction_prob, copy=False)
        death_prob = 1 / (1 + numpy.exp(
            fitness_grid * death_parameter_one - death_parameter_two))
        numpy.nan_to_num(death_prob, copy=False)

        population_grid += numpy.random.binomial(population_grid, reproduction_prob)
        population_grid -= numpy.random.binomial(population_grid, death_prob)

        decisions = [
            i
            for i, n_decisions in numpy.ndenumerate(numpy.random.binomial(
                    population_grid, decision_making_probability))
            for _ in range(n_decisions)]

        while decisions:
            cell = decisions.pop(numpy.random.randint(len(decisions)))
            fitness = fitness_grid[cell]

            copy_candidates = []
            neighbors = list(neighborhood(cell))
            for n in neighbors:
                if numpy.random.binomial(population_grid[n], observed_agents):
                    copy_candidates.append((fitness_grid[n], numpy.random.random(), n))
            target_fitness, _, target_cell = max(
                copy_candidates,
                default=(None, None, None))
            pop = population_grid[cell]
            target_pop = 0 if target_cell is None else population_grid[target_cell]

            # Follow the more concise description of Crema (2015) in the decision process
            if pop > 1 and target_pop > 1:
                if fitness <= basic_individual_payoff - evidence_threshold and (
                        target_fitness <= basic_individual_payoff - evidence_threshold or
                        fitness >= target_fitness):
                    fission(cell, neighbors, population_grid)
                elif target_fitness > basic_individual_payoff - evidence_threshold and (
                        fitness <= target_fitness - evidence_threshold or
                        fitness <= basic_individual_payoff - evidence_threshold):
                    population_grid[cell] -= 1
                    population_grid[target_cell] += 1
            elif pop > 1 and target_pop == 1:
                if fitness < basic_individual_payoff - evidence_threshold or fitness < target_fitness - evidence_threshold:
                    fission(cell, neighbors, population_grid)
            elif pop == 1 and target_pop > 1:
                if fitness < target_fitness - evidence_threshold:
                    population_grid[cell] -= 1
                    population_grid[target_cell] += 1
            elif pop == 1 and target_pop == 1:
                if fitness < basic_individual_payoff > target_fitness:
                    # The description contains a typo here, stating 'Fission' but meaning 'Fusion'
                    population_grid[cell] -= 1
                    population_grid[target_cell] += 1
                    try:
                        decisions.remove(target_cell)
                    except ValueError:
                        pass
            elif pop > 1 and target_cell is None:
                if fitness < basic_individual_payoff - evidence_threshold:
                    fission(cell, neighbors, population_grid)
        run_results.append(population_grid + 0)

    return run_results

if __name__ == '__main__':
    from drennan2004comparing import a_coefficient
    import csv
    results = {
        "A coefficient mean over the last 200 steps": lambda a: numpy.mean(a[-200:]),
        "A coefficient std over the last 200 steps": lambda a: numpy.std(a[-200:]),
        "Maximum A coefficient over the last 200 steps": lambda a: max(a[-200:]),
        "Minimum A coefficient over the last 200 steps": lambda a: min(a[-200:]),
    }
    with open("summary.csv", "w") as csvfile:
        csvoutput = csv.DictWriter(
                csvfile,
                list(itertools.chain(["run"], model_parameters.keys(), results.keys())))
        csvoutput.writeheader()
        for i, parameters in enumerate(itertools.product(*model_parameters.values())):
            for j in range(10):
                run_id = "run_{:05d}{:01d}".format(i, j)
                named_parameters = dict(zip(model_parameters.keys(), parameters))
                named_parameters["neighborhood"] = lambda x: neighbors_within(
                    x, named_parameters["spatial_interaction_range"] * 1.5, (10, 10), (True, True))
                print(named_parameters)
                run_results = simulation(**named_parameters)
                del named_parameters["neighborhood"]
                json.dump(named_parameters, open("{:s}-parameters.json".format(run_id), "w"),
                        indent=2,
                        sort_keys=True)
                a_coefficients = [a_coefficient([pop for pop in x.flat if pop])
                                for x in run_results]
                plt.plot(a_coefficients)
                named_parameters.update({
                    key: val(a_coefficients)
                    for key, val in results.items()})
                named_parameters["run"] = run_id
                csvoutput.writerow(named_parameters)
                csvfile.flush()
                plt.ylim(-1.6, 1.1)
                plt.savefig("{:s}.png".format(run_id))
                plt.close()
