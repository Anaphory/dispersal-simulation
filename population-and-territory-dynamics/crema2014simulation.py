"""crema2014simulation

A replication of the model described by Crema (2014).

Crema, Enrico R. 2014. A Simulation Model of Fission–Fusion Dynamics and
Long-Term Settlement Change. Journal of Archaeological Method and Theory 21(2).
385–404. (doi:10.1007/s10816-013-9185-4)

"""
import numpy
from matplotlib import pyplot as plt

model_parameters = dict(
    initial_agents = 10,
    time_steps = 500,
    patches = 100,
    basic_individual_payoff = 10,
    payoff_standard_deviation = 1.,
    # The paper calls this parameter ε “payoff variance”, but also describes
    # the process as “random draw from a Gaussian probability distribution with
    # […] a standard deviation ε”. Because the value is 1, it does not actually
    # make any difference.
    cooperative_benefit = [0.3, 0.5, 0.8],
    resource_input = 200,
    basic_reproduction_rate = 0.05,
    death_parameter_one = [0.8, 1.0, 1.2, 1.4],
    death_parameter_two = 5,
    spatial_interaction_range = [1, 100],
    decision_making_probability = [0.1, 0.5, 1.0],
    evidence_threshold = 3,
    observed_agents = [1e-7, 0.5, 1])


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
        The neigboring cells, from which an empty one is drawn at random
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
        print("Attempt to fission from {:} failed".format(cell))
        return
    modifyable_population_grid[cell] -= 1
    target = empty_cells[numpy.random.randint(len(empty_cells))]
    assert population_grid[target] == 0
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
        whether to wrap the neigborhood around the borders. (False, False)
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
        neigborhood = lambda x: neighbors_within(
            x, 1.5, (10, 10), (True, True)),
        spatial_interaction_range = 1.5,
        decision_making_probability = 0.1,
        evidence_threshold = 3,
        observed_agents = 1e-7):
    """Run Crema's simulation

    Returns
    =======
    list(numpy.array)
        The list of the state of the population in each time step
    """
    grid_shape = int(patches ** 0.5), patches // int(patches ** 0.5)
    population_grid = numpy.zeros(grid_shape, dtype=int)
    resource_grid = resource_input * numpy.ones_like(population_grid)
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
                foraging_contribution, resource_grid) / population_grid

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
            for i, n_decisions in enumerate(numpy.random.binomial(
                    population_grid, decision_making_probability).flat)
            for _ in range(n_decisions)]
        while decisions:
            i = decisions.pop(numpy.random.randint(len(decisions)))
            cell = numpy.unravel_index(i, population_grid.shape)
            copy_candidates = []
            neighbors = list(neigborhood(cell))
            for n in neighbors:
                if numpy.random.binomial(population_grid[n], observed_agents):
                    copy_candidates.append((fitness_grid[n], numpy.random.random(), n))
            other_fitness, _, target_cell = max(
                copy_candidates,
                default=(None, None, None))
            if population_grid[cell] == 1:
                if target_cell is None:
                    pass
                elif  population_grid[target_cell] == 1:
                    if fitness_grid[cell] >= other_fitness:
                        pass
                    elif (fitness_grid[cell] < basic_individual_payoff and
                        other_fitness < basic_individual_payoff):
                        population_grid[cell] -= 1
                        population_grid[target_cell] += 1
                else:
                    if fitness_grid[cell] <= other_fitness - evidence_threshold:
                        population_grid[cell] -= 1
                        population_grid[target_cell] += 1
            else:
                if target_cell is None:
                    if fitness_grid[cell] <= basic_individual_payoff - evidence_threshold:
                        fission(cell, neighbors, population_grid)
                elif population_grid[target_cell] > 1:
                    if other_fitness > basic_individual_payoff - evidence_threshold:
                        if fitness_grid[cell] <= max(basic_individual_payoff, other_fitness):
                            population_grid[cell] -= 1
                            population_grid[target_cell] += 1
                    elif other_fitness <= basic_individual_payoff - evidence_threshold:
                        if fitness_grid[cell] <= basic_individual_payoff - evidence_threshold:
                            fission(cell, neighbors, population_grid)
                    elif fitness_grid[cell] >= other_fitness:
                        if fitness_grid[cell] > basic_individual_payoff - evidence_threshold:
                            pass
                        else:
                            fission(cell, neighbors, population_grid)
                else: # population_grid[target_cell] == 1:
                    if fitness_grid[cell] >= other_fitness:
                        if fitness_grid[cell] > basic_individual_payoff - evidence_threshold:
                            pass
                        else:
                            fission(cell, neighbors, population_grid)
                    elif fitness_grid[cell] <= max(
                            other_fitness, basic_individual_payoff) - evidence_threshold:
                        fission(cell, neighbors, population_grid)

        run_results.append(population_grid + 0)

    return run_results

if __name__ == '__main__':
    run_results = simulation()
    plt.plot([x.sum() for x in run_results])
    plt.show()
