"""delcastillo2013modeling

Implement a geographic (abstracted) agent_based model of migrating
hunter_gatherer bands with culture. The implementation is based on the NetLogo
implementation of PSMED, published by Baceló et al. (2013), while making sure
to implement the specifics of del Castillo (2013).

Barceló, Joan A. & Cuesta, J A & del Castillo, Florencia & del Olmo, Ricardo &
    Galán, José M & Mameli, Laura & Miguel, Francisco J & Poza, David & Santos,
    José I & Vilà, Xavier. 2013. PSMED - Patagonia Simple Model of Ethnic
    Differentiation. CoMSES Computational Model Library. (4063/releases/1.0.0,
    CoMSES.) (https://www.comses.net/codebases/4063/releases/1.0.0/) (Accessed
    February 17, 2020.)

del Castillo, F. & Barceló, J. A. & Mameli, L. & Miguel, F. & Vila, X. 2013.
    Modeling Mechanisms of Cultural Diversity and Ethnicity in Hunter–Gatherers.
    Journal of Archaeological Method and Theory 21(2). 364–384.
    (doi:10.1007/s10816-013-9199-y)

"""

import attr
import numpy
import networkx as nx
import collections
import matplotlib.pyplot as plt

@attr.s
class Family:
    collected_energy = attr.ib()      # energy obtained during the current tick
    identity = attr.ib(eq=False)              # 10-dimension vector with values 1 (no important) to 6 (very important)
    labor = attr.ib()                 # number of individuals within a household
    technology = attr.ib()            # Gaussian_distributed: mean: average_technology; standard deviation: diversity

    my_neighborhood = attr.ib(factory=list, repr=False)       # list of families around
    my_helpers = attr.ib(factory=list, repr=False)            # subof = my_neighborhood that will help me (self,with probability 95%)
    my_group = attr.ib(factory=list, repr=False)              # list of agents withins a social aggregate (self,i.e. families that belong to my group/component)
    group_size = attr.ib(default=0)            # size of the group (self,component) where I belong
    total_energy = attr.ib(default=0)          # total level of energy (self,including surplus)
    individual_capability = attr.ib(default=1) # fraction of the resource that the agent is able to obtain when working on their own
    cultural_distance = attr.ib(default=0)     # cooperation will be posible if similarity >= cultural_distance
    cooperation = attr.ib(default=False)           # True if I have been involved in a cooperation process (self,i.e. I have helped someone or someone has helped me)
    explored = attr.ib(default=False)              # needed for detecting components
    component = attr.ib(default=0)             # group (self,i.e. component in the network) where the family belongs
    _patch = attr.ib(default=None, repr=False)
    _links = attr.ib(factory=list, repr=False)

    @property
    def conservation_factor(self):
        # current value of technology divided by two
        return self.technology / 2

    @property
    def survival_threshold(self):
        # 730 * labor
        return (730 * self.labor) / 2

    @property
    def patch(self):
        return self._patch

    @patch.setter
    def patch(self, patch):
        if self._patch is not None:
            self._patch.families.remove(self)
        self._patch = patch
        self._patch.families.append(self)

    def linked_neighbors(self):
        return self._links

    def create_link_with(self, other):
        if other == self:
            return
        self._links.append(other)
        other._links.append(self)

    def copy_to(self, patch):
        new = Family(
            collected_energy = 0,
            identity = self.identity + 0,
            labor = self.labor,
            technology = self.technology,
        )
        new.patch = patch
        new.create_link_with(self)

@attr.s
class Patch:
    difficulty = attr.ib()      # (h)
    max_resource = attr.ib()    # maximum amount of resource (in hot season)
    resource = attr.ib()        # current amount of resource
    exploited = attr.ib()       # True during one tick after been exploited
    coords = attr.ib()
    families = attr.ib(factory=list, repr=False)


def distance(patch1, patch2):
    return ((patch1.coords[0] - patch2.coords[0]) ** 2 + (patch1.coords[1] - patch2.coords[1]) ** 2) ** 0.5


class Simulation():
    def __init__(self):
        self.steps = 10000 # 100
        self.diversity = 0.5
        self.labor_average = 2
        self.internal_change_rate = 0.05
        self.initial_population = 20 # 100
        self.average_technology = 0.22
        self.movement = 5
        self.max_resource_on_patches = 20000
        self.cooperation_allowed = True

        self.number_of_agents_that_died_of_starvation = 0

    def setup(self):
        #INITIALIZE PATCHES
        self.patches = [
            Patch(difficulty = numpy.random.random(),
                  # Resources are uniformly distributed in [max_resource_on_patches / 1000, max_resource_on_patches]
                  max_resource = (self.max_resource_on_patches / 1000) + numpy.random.random() * ( self.max_resource_on_patches - (self.max_resource_on_patches / 1000) ),
                  resource = 10,
                  exploited = False,
                  coords=(i, j))
            for i in range(50) for j in range(50)]
        self.tick = 0

        #GLOBAL = VARIABLES
        self.weights_vector = numpy.array([0.5, 0.25, 0.12, 0.06, 0.03, 0.015, 0.0075, 0.00375, 0.001875, 0.0009375])

        #CREATE AND INITIALIZE FAMILIES
        initial_identity = numpy.random.randint(1, 7, size=len(self.weights_vector))  #all the families start with the same identity vector

        self.families = []
        randomstate = numpy.random.default_rng(0)
        for patch in randomstate.choice(self.patches, self.initial_population, replace=False):
            labor = numpy.random.poisson(self.labor_average)
            if labor < 2:
                labor = 2 #(correction for cases with labor < 2)
            family = Family(
                identity = initial_identity,
                labor = labor,
                technology = max(min(numpy.random.normal(self.average_technology, self.diversity), 2), 0.02),
                        # correction to force technology to be in the interval [0 2] (and thus depreciation will be in [0 1]
                collected_energy = 0,
                my_group = []     # initialize my_group = nobody
                )
            family.total_energy = 1 + family.survival_threshold * (1 + family.conservation_factor) / family.conservation_factor
            family.patch = patch
            # this means that all the families will have enough energy to survive the first tick without hunting
            self.families.append(family)


    def survive(self):
        self.number_of_living_agents = len(self.families)

        # 1. Substract survival_threshold. If they didn't reach their survival
        # threshold in the previous tick, their current level of total_energy will be
        # < 0 at this point; in this case their amount of labor is reduced by one
        # unit
        for family in self.families:
            family.total_energy -= family.survival_threshold

            if family.total_energy < 0:
                family.labor -= 1
                family.total_energy = 0

                #if their amount of labor goes below 2, the family (agent) dies
                if family.labor < 2:
                    self.number_of_agents_that_died_of_starvation += 1
                    self.families.remove(family)
                    family.patch.families.remove(family)

        # 2. Vegatative reproduction: The amount of labor within a family goes up one
        # unit every 6 ticks – The original code actually says ‘30’ instead of ‘6’
        if self.tick % 6 == 0:
            for family in self.families:
                family.labor += 1

        # if the amount of labor reaches the value 10, the family is split into two
        # families with probability 95% (provided there is enough room in the
        # neighbourhood for another family: there cannot be more than one family in
        # the same patch)
        for family in self.families:
            if family.labor >= 10:
                if numpy.random.random() < .95:
                    reachable_empty_patches = list(self.reachable_empty_patches_within(family.patch, self.movement))
                    if reachable_empty_patches:
                            family.labor //= 2
                            family.copy_to(reachable_empty_patches[numpy.random.randint(len(reachable_empty_patches))])
    def reachable_empty_patches_within(self, patch, radius):
        for other in self.patches:
            if other.families:
                continue
            if distance(patch, other) < radius:
                yield other

    def update_resources_on_patches(self):
        season = ["hot", "cold"][self.tick % 2]
        for patch in self.patches:
            if patch.exploited:
                # patches that were exploited in the previous tick:
                # if the current season is cold, the patch will not regenerate completely: (max_resource - amount_consumed)/2
                # if the current season is hot, the level of resource will be max_resource (no matter whether it was exploited during the previous tick or not)
                if season == "hot":
                    patch.resource = patch.max_resource
                else:
                    patch.resource = (patch.max_resource - patch.amount_consumed) / 2
            else:
                patch.resource = patch.max_resource if (season == "hot") else (patch.max_resource / 2)
            patch.exploited = False
            patch.amount_consumed = 0


    def update_cultural_distance(self):
        for family in self.families:
            if (family.patch.resource) <= family.survival_threshold:
                # They will refuse to cooperate, as they are unable to reach their survival_threshold no matter the number of people who help them
                family.cultural_distance = 1
            else:
                extra_workers = max(
                    0,
                    ((family.survival_threshold / (family.patch.difficulty * (family.patch.resource - family.survival_threshold))) ** (1 / family.technology) - family.labor))
                # extra_workers could be negative (if I have plenty of labor to obtain my survival threshold: I don't need any help from anyone)  --> in this case extra_workers = 0
                family.cultural_distance = min(
                    1, ((1 / 100) * extra_workers)) # cultural distance could be greater than 1 if the number of extra workers goes above 100 --> --> in this case cultural distance = 1


    def hunt_and_gather(self):
        for family in self.families:
            family.total_energy = family.total_energy * family.conservation_factor #depreciate agents' energy at the beginning of the tick
            family.cooperation = False
            family.my_helpers = [] # (empty agentlist)
            family.my_group = []   # (empty agentlist)
            family.collected_energy = 0

        # families with total_energy > survival_threshold don't hunt or gather

        for family in self.families:
            if not family.total_energy < family.survival_threshold:
                continue
            # Try to act individually:
            family.individual_capability = self.calculate_individual_capability(family)
            family.productivity = ( family.patch.resource) * family.individual_capability

            if self.cooperation_allowed:
                # FAMILIES ARE ALLOWED TO ASK FOR HELP WHENEVER THEY ARE UNABLE TO REACH THEIR SURVIVAL THRESHOLD
                if family.productivity >= family.survival_threshold:
                    # ACT INDIVIDUALLY
                    # I don't need to cooperate: I will act individually and collect as much energy as I need to reach my survival_threshold
                    family.collected_energy = min((family.survival_threshold - family.total_energy), (family.individual_capability) * family.patch.resource)
                    family.total_energy = family.total_energy + family.collected_energy
                    p = family.patch #update resource on patch:
                    amount_consumed = family.collected_energy
                    p.resource = p.resource - amount_consumed
                    p.exploited = True
                else:
                    # COOPERATION (I need to ask for help)
                    family.capability = family.individual_capability
                    self.identify_neighbors(family) # define my_helpers (self,my_helpers contains a list of families that can potentially help me, i.e. our similarity is greater than their cultural_distance)
                    agents_willing_to_help = self.ask_for_help(family)
                    if agents_willing_to_help: # if someone helps me, my capability will be aggregated_capability. Otherwise it will be my individual_capability
                        family.capability = self.calculate_aggregated_capability(family, agents_willing_to_help)

                    family.collected_energy = (family.patch.resource) * family.capability

                    family.total_energy = family.total_energy + family.collected_energy        #(therefore, total_energy might be greater than her survival_threshold

                    #update resource on patch:
                    p = family.patch
                    amount_consumed = family.collected_energy
                    p.resource = p.resource - amount_consumed
                    p.exploited = True
                # FAMILIES ARE NOT ALLOWED TO ASK FOR HELP IF THEY ARE UNABLE TO REACH THEIR SURVIVAL THRESHOLD. THEY WILL ACT INDIVIDUALLY
                family.collected_energy = min((family.survival_threshold - family.total_energy), (family.individual_capability) * family.patch.resource)
                family.total_energy = family.total_energy + family.collected_energy
                p = family.patch #update resource on patch:
                amount_consumed = family.collected_energy
                p.resource = p.resource - amount_consumed
                p.exploited = True


    def calculate_individual_capability(self, family):
        return 1 / (1 + (1 / ((family.patch.difficulty) * (family.labor ** family.technology))))


    def ask_for_help(self, myself):
        # Determine the list of agents that are willing to cooperate
        agents_willing_to_cooperate = []
        if myself.my_helpers:  # if my list of helpers is not empty
            #each one of my helpers will help me with p=95%
            for helper in myself.my_helpers:
                if numpy.random.random() < 0.95:
                    agents_willing_to_cooperate.append(helper)
            # The agents contained in the list agents_willing_to_cooperate will help me
            for agent in agents_willing_to_cooperate:
                print("Creating link:", myself.patch.coords, agent.patch.coords, distance(agent.patch, myself.patch))
                agent.create_link_with(myself)   #helpers create a link with helped agents
            myself.cooperation = True # I have cooperated...
            for agent in agents_willing_to_cooperate:
                agent.cooperation = True  # ... and so have my helpers
        return agents_willing_to_cooperate


    def calculate_aggregated_capability(self, family, agents_willing_to_cooperate):
        max_technology_within_my_local_group = max(family.technology for family in agents_willing_to_cooperate)   # maximum technology within my group (including myself)
        returns_to_cooperation = 2 - (family.patch.resource / self.max_resource_on_patches)
        total_labor = family.labor + sum(family.labor for family in agents_willing_to_cooperate) # total labor within my group (including my)self)
        aggregated_capability = 1 / (1 + ( 1 / ( ( (family.patch.difficulty) * (total_labor ** max_technology_within_my_local_group)) ** returns_to_cooperation ) ) )
        return aggregated_capability


    def identify_neighbors(self, myself):
        myself.my_neighborhood = [family
                                  for family in self.families
                                  if distance(family.patch, myself.patch) <= self.movement]
        print(myself.patch.coords, [family.patch.coords for family in myself.my_neighborhood])
        print([distance(myself.patch, family.patch) for family in myself.my_neighborhood])
        myself.my_helpers = [family
                      for family in myself.my_neighborhood
                      if self.get_similarity(family.identity, myself.identity) > myself.cultural_distance] # cultural_distance of my neighbors'
        # FIXME: YES, THE ORIGINAL SOURCE HAS A > THERE


    def get_similarity(self,a, b):
        ap = a * self.weights_vector
        bp = b * self.weights_vector
        numerator = (ap * bp).sum()
        denominator = (ap * ap).sum() ** 0.5 * (bp * bp).sum() ** 0.5
        return numerator / denominator


    def identify_groups(self):
        self.find_all_components()
        for myself in self.families:
            if not myself.component > 0:
                continue
            myself.my_group = [family for family in self.families if family.component == myself.component]
            myself.group_size = len(myself.my_group)


    ## The following two procedures are an adaptation from: Wilensky, U. (2005). NetLogo Giant Component model.
    ## http://ccl.northwestern.edu/netlogo/models/GiantComponent.
    ## Center for Connected Learning and Computer_Based Modeling, Northwestern University, Evanston, IL.
    def find_all_components(self):
        for myself in self.families:
            myself.explored = False
            myself.group_size = 0
        for myself in self.families:
            if not myself.linked_neighbors():
                myself.component = 0
                myself.explored = True # families that don't cooperate (isolated agents) will have component = 0
        self.component_index = 0
        while True:
            if [family for family in self.families if not family.explored]:
                    start = numpy.random.choice([family for family in self.families if not family.explored])
                    self.component_index = ( self.component_index + 1)
                    self.explore(start, self.component_index)
            else:
                break


    ## Finds all families reachable from this node
    def explore(self, myself, component_index):
        if myself.explored:
            return
        myself.explored = True
        myself.component = self.component_index
        for neighbor in myself.linked_neighbors():
            self.explore(neighbor, component_index)


    def decide_whether_to_move_or_stay(self):
        for myself in self.families:
            # 1. Calculate if I will be able to survive one more tick with my current level of energy (bearing in mind the depreciation factor) without working.
            # If that is the case, I leave with probability 0.05
            if (myself.total_energy - myself.survival_threshold) * myself.conservation_factor > myself.survival_threshold:
                if (numpy.random.random() < 0.05):
                    move()
                # 2. I will have to work next tick because I will not be able to get by with my (depreciated) level of energy
                # Before moving, I will check if I will be able to get enough resources if I stay in the same patch
                myself.patch.resources_next_tick =  (myself.patch.max_resource - collected_energy) / 2 if (season == "hot") else myself.patch.max_resource
                if resources_next_tick * individual_capability  > survival_threshold:
                    # If I can survive here, I will probably stay, but I will leave with probability 0.05
                    if (numpy.random.random() < 0.05):
                        move()
                else:
                    move()


    def move(self):
        if [patch for patch in patches if distance(patch, this_family.patch) <= movement if not patch.families]:
            self.number_of_movements = self.number_of_movements + 1
            move_to(numpy.random.choice([patch for patch in patches if distance(patch, this_family.patch) <= movement if not patch.families]))


    def update_identity(self):
        # 1. Diffusion process (only for agents that have cooperated)
        if self.families:
            group_index = range(1, self.component_index + 1)
            for g in group_index:
                consensual_identity = self.compute_consensual_identity([family for family in self.families if family.component == g])
                for family in [family for family in self.families if family.component == g]:
                    if numpy.random.random() < 0.95:
                        family.identity = consensual_identity
        # 2. Mutation process (for all families)
        for myself in self.families:
            if numpy.random.random() < self.internal_change_rate:
                index_vector = range(10)
                for i in index_vector:
                        if numpy.random.random() >= (self.weights_vector[i]):
                            myself.identity[i] = numpy.random.randint(1, 7)


    def compute_consensual_identity(self,group):
        consensus = numpy.zeros(10, int)
        for n in range(10):
                consensual_trait = numpy.random.choice([family.identity[n] for family in group])
                consensus[n] = consensual_trait
        return consensus


    def update_technology(self):
        if self.families:
            group_index = range(1, self.component_index + 1)
            for g in group_index:
                this_group = [family for family in self.families if family.component == g]
                average_technology_this_group = numpy.mean([family.technology for family in this_group])
                for myself in this_group:
                    if myself.technology < 0.95 * average_technology_this_group:
                        if numpy.random.random() < 0.95:
                            myself.technology = myself.technology + 0.1
                            if myself.technology > max([family.technology for family in this_group]):
                                myself.technology = max([family.technology for family in this_group])
                        if myself.technology < 1.05 * average_technology_this_group:
                            if numpy.random.random() < 0.95:
                                myself.technology = myself.technology + 0.01
                                if myself.technology > max([family.technology for family in this_group]):
                                    myself.technology = max([family.technology for family in this_group])

        # 2. Mutation process (for all families)
        for myself in self.families:
            if numpy.random.random() < self.internal_change_rate:
                myself.technology = numpy.random.normal(self.average_technology, self.diversity)
                # possible correction to force technology to be in the interval [0 2] (and thus depreciation will be in [0 1]
                if myself.technology < 0.02:
                    myself.technology = 0.02
                if myself.technology > 2:
                    myself.technology = 2


    def update_output_variables(self):
        number_of_social_aggregates = self.component_index
        number_of_agents_in_social_aggregates = len([family for family in self.families if family.cooperation])
        number_of_agents_in_social_aggregates = len([family for family in self.families if not family.cooperation])
        largest_group_size = max([family.group_size for family in self.families], default=None)
        if largest_group_size == 1:
            largest_group_size = "N/A"
        total_collected_energy = sum([family.collected_energy for family in self.families])
        std_collected_energy = numpy.std([family.collected_energy for family in self.families])
        average_cultural_distance_in_aggregates = numpy.mean(
            [family.cultural_distance for family in self.families if family.component > 0])
        average_cultural_distance_in_aggregates = numpy.std(
            [family.cultural_distance for family in self.families if family.component > 0])
        sum_of_labor = sum([family.labor for family in self.families])
        if len( self.families) >= 2:
            mean_technology_of_families = numpy.mean([family.technology for family in self.families])
            std_technology_of_families = numpy.std([family.technology for family in self.families])
        print(locals())

    def plot(self):
        resources = numpy.zeros((50, 50))
        for p in self.patches:
            resources[p.coords] = p.resource
        plt.imshow(resources)
        g = nx.Graph()
        for family in self.families:
            g.add_node(family.patch.coords)
            for other in family.linked_neighbors():
                g.add_edge(family.patch.coords, other.patch.coords)
        nx.draw_networkx_nodes(g, mirror_dict(),
                               node_color = "#dd0000",
                               node_size = 10,
                               node_shape = "d")
        nx.draw_networkx_edges(g, mirror_dict())
        plt.savefig("tick-{:09d}.png".format(self.tick))
        plt.close()

    def simulate(self):
        for step in range(self.steps):
            self.survive()
            self.update_resources_on_patches()
            self.update_cultural_distance()
            self.hunt_and_gather()
            self.identify_groups()
            self.update_identity()
            self.update_technology()
            self.decide_whether_to_move_or_stay()
            self.update_output_variables()
            self.plot()
            self.tick += 1
            if not self.families:
                break


class mirror_dict(collections.defaultdict):
    def __missing__(self, key):
        return key

for i in range(20):
    s = Simulation()
    s.setup()
    s.simulate()
