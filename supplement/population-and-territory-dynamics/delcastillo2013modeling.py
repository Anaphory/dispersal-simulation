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

if False:
  season                         # hot or cold
  weights_vector                 # [0.5 0.25 0.12 0.06 0.03 0.015 0.0075 0.00375 0.001875 0.0009375]
  component_index                # counter of groups

  #OUTPUT VARIABLES
  number_of_social_aggregates
  number_of_agents_in_social_aggregates
  number_of_isolated_agents
  number_of_agents_that_died_of_starvation
  number_of_living_agents      #(at the beginning of the tick)
  total_collected_energy
  collected_energy_standard_deviation
  average_cultural_distance_in_aggregates
  sd_cultural_distance_in_aggregates
  largest_group_size
  number_of_movements
  total_number_of_starvation_deaths
  sum_of_labor
  mean_technology_of_families
  std_technology_of_families


@attr.s
class Family:
  my_neighborhood = attr.ib()       # list of families around
  my_helpers = attr.ib()            # subof = my_neighborhood that will help me (with probability 95%)
  my_group = attr.ib()              # list of agents withins a social aggregate (i.e. families that belong to my group/component)
  group_size = attr.ib()            # size of the group (component) where I belong
  collected_energy = attr.ib()      # energy obtained during the current tick
  total_energy = attr.ib()          # total level of energy (including surplus)
  identity = attr.ib()              # 10-dimension vector with values 1 (no important) to 6 (very important)
  labor = attr.ib()                 # number of individuals within a household
  survival_threshold = attr.ib()    # 730 * labor
  technology = attr.ib()            # Gaussian_distributed: mean: average_technology; standard deviation: diversity
  conservation_factor = attr.ib()   # current value of technology divided by two
  individual_capability = attr.ib() # fraction of the resource that the agent is able to obtain when working on their own
  cultural_distance = attr.ib()     # cooperation will be posible if similarity >= cultural_distance
  cooperation = attr.ib()           # true if I have been involved in a cooperation process (i.e. I have helped someone or someone has helped me)
  explored = attr.ib()              # needed for detecting components
  component = attr.ib()             # group (i.e. component in the network) where the family belongs


@attr.s
class Patch:
  difficulty = attr.ib()      # (h)
  max_resource = attr.ib()    # maximum amount of resource (in hot season)
  resource = attr.ib()        # current amount of resource
  amount_consumed = attr.ib() # amount of energy consumed during the current tick
  exploited = attr.ib()       # true during one tick after been exploited


def setup():
  clear_all()
  reset_ticks()

  #GLOBAL = VARIABLES
  weights_vector = [0.5, 0.25, 0.12, 0.06, 0.03, 0.015, 0.0075, 0.00375, 0.001875, 0.0009375]

  #CREATE AND INITIALIZE FAMILIES
  initial_identity = numpy.random.randint(1, 7)  #all the families start with the same identity vector

  randomstate = np.random.default_rng()
  for i in randomstate.choice(patches, initial_population, replace=False):
    labor = numpy.random.poisson(labor_average)
    if labor < 2:
      labor = 2 #(correction for cases with labor < 2)
    Family(
      shape = "person",
      color = red,
      identity = initial_identity,
      labor = labor,
      survival_threshold = (730 * labor) / 2,
      technology = max(min(random_normal( average_technology, diversity), 2), 0.02),
          # correction to force technology to be in the interval [0 2] (and thus depreciation will be in [0 1]
      conservation_factor = technology / 2,
      total_energy = 1 + survival_threshold * (1 + conservation_factor) / conservation_factor,    # this means that all the families will have enough energy to survive the first tick without hunting
      collected_energy = 0,
      my_group = []     # initialize my_group = nobody
      )

  #INITIALIZE PATCHES
  for patch in patches:
    # Resources are uniformly distributed in [max_resource_on_patches / 1000, max_resource_on_patches]
    max_resource = (max_resource_on_patches / 1000) + random ( max_resource_on_patches - (max_resource_on_patches / 1000) )
    difficulty = numpy.random.random()
    exploited = false
    amount_consumed = 0
    pcolor = 69

def survive(families):
  number_of_agents_that_died_of_starvation = 0
  number_of_living_agents = len(families)

  # 1. Substract survival_threshold. If they didn't reach their survival
  # threshold in the previous tick, their current level of total_energy will be
  # < 0 at this point; in this case their amount of labor is reduced by one
  # unit
  for family in families:
      family.total_energy -= family.survival_threshold

      if family.total_energy < 0:
        family.labor -= 1
        family.survival_threshold = (730 * family.labor) / 2
        family.total_energy = 0

        #if their amount of labor goes below 2, the family (agent) dies
        if family.labor < 2:
          number_of_agents_that_died_of_starvation += 1
          family.die()

  # 2. Vegatative reproduction: The amount of labor within a family goes up one
  # unit every 6 ticks – The code actually says ‘30’ instead of ‘6’
  for family in families:
    if ticks % 30 == 0:
      family.labor += 1
      family.survival_threshold = (730 * family.labor) / 2

  # if the amount of labor reaches the value 10, the family is split into two
  # families with probability 95% (provided there is enough room in the
  # neighbourhood for another family: there cannot be more than one family in
  # the same patch)
  for family in families:
    if family.labor >= 10:
      if numpy.random.random() < .95:
        reachable_empty_patches = reachable_empty_patches_within(family.location, in_radius)
        if reachable_empty_patches:
            family.labor //= 2
            family.survival_threshold = (730 * family.labor) / 2
            family.copy_to(reachable_empty_patches[numpy.random.randint(len(reachable_empty_patches))])


def update_resources_on_patches():
  season = ["hot", "cold"][ticks % 2]
  for patch in patches:
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


def update_cultural_distance():
  for family in families:
    if (family.patch.resource) <= family.survival_threshold:
      # They will refuse to cooperate, as they are unable to reach their survival_threshold no matter the number of people who help them
      family.cultural_distance = 1
    else:
      extra_workers = max(
        0,
        ((survival_threshold / (difficulty * (family.patch.resource - survival_threshold))) ** (1 / technology) - labor))
      # extra_workers could be negative (if I have plenty of labor to obtain my survival threshold: I don't need any help from anyone)  --> in this case extra_workers = 0
      cultural_distance = min(
        1, ((1 / 100) * extra_workers)) # cultural distance could be greater than 1 if the number of extra workers goes above 100 --> --> in this case cultural distance = 1


def hunt_and_gather():
  for family in families:
    total_energy = total_energy * conservation_factor #depreciate agents' energy at the beginning of the tick
    cooperation = false
    my_helpers = [] # (empty agentlist)
    my_group = []   # (empty agentlist)
    collected_energy = 0

  links = []


  # families with total_energy > survival_threshold don't hunt or gather

  for family in families:
    if not total_energy < survival_threshold:
      continue
    # Try to act individually:
    individual_capability = calculate_individual_capability
    productivity = ( family.patch.resource) * individual_capability

    if cooperation_allowed:
      # FAMILIES ARE ALLOWED TO ASK FOR HELP WHENEVER THEY ARE UNABLE TO REACH THEIR SURVIVAL THRESHOLD
      if productivity >= survival_threshold:
        # ACT INDIVIDUALLY
        # I don't need to cooperate: I will act individually and collect as much energy as I need to reach my survival_threshold
        collected_energy = min ( (survival_threshold - total_energy), (individual_capability) * family.patch.resource)
        total_energy = total_energy + collected_energy
        with family.patch: #update resource on patch:
            amount_consumed = family.collected_energy
            resource = resource - amount_consumed
            exploited = true
      else:
        # COOPERATION (I need to ask for help)
        capability = individual_capability
        identify_neighbors # define my_helpers (my_helpers contains a list of families that can potentially help me, i.e. our similarity is greater than their cultural_distance)
        agents_willing_to_help = ask_for_help
        if agents_willing_to_help: # if someone helps me, my capability will be aggregated_capability. Otherwise it will be my individual_capability
          capability = calculate_aggregated_capability(agents_willing_to_help)

        collected_energy = (family.patch.resource) * capability

        total_energy = total_energy + collected_energy        #(therefore, total_energy might be greater than her survival_threshold

        #update resource on patch:
        with family.patch:
            amount_consumed = family.collected_energy
            resource = resource - amount_consumed
            exploited = true
      # FAMILIES ARE NOT ALLOWED TO ASK FOR HELP IF THEY ARE UNABLE TO REACH THEIR SURVIVAL THRESHOLD. THEY WILL ACT INDIVIDUALLY
      collected_energy = min ( (survival_threshold - total_energy) (individual_capability) * family.patch.resource)
      total_energy = total_energy + collected_energy
      with family.patch: #update resource on patch:
            amount_consumed = family.collected_energ
            resource = resource - amount_consumed
            exploited = true


def calculate_individual_capability():
  return 1 / (1 + (1 / ((family.patch.difficulty) * (labor ** technology))))


def ask_for_help():
  # Determine the list of agents that are willing to cooperate
  agents_willing_to_cooperate = []
  if my_helpers:  # if my list of helpers is not empty
    #each one of my helpers will help me with p=95%
    for helper in my_helpers:
      if random() < 0.95:
        agents_willing_to_cooperate.append(helper)
    # The agents contained in the list agents_willing_to_cooperate will help me
    for agent in turtle_agents_willing_to_cooperate:
      agent.create_link_with(family)   #helpers create a link with helped agents
    cooperation = true # I have cooperated...
    for agent in turtle_agents_willing_to_cooperate():
      agent.cooperation = true  # ... and so have my helpers
  return agents_willing_to_cooperate


def calculate_aggregated_capability(agents_willing_to_cooperate):
  max_technology_within_my_local_group = max(family.technology for family in agents_willing_to_cooperate())   # maximum technology within my group (including myself)
  returns_to_cooperation = 2 - (family.patch.resource / max_resource_on_patches)
  total_labor = family.labor + sum(family.labor for family in agents_willing_to_cooperate()) # total labor within my group (including myself)
  aggregated_capability = 1 / (1 + ( 1 / ( ( (family.patch.difficulty) * (total_labor ** max_technology_within_my_local_group)) ** returns_to_cooperation ) ) )
  return aggregated_capability


def identify_neighbors():
  my_neighborhood = [family for family in families if distance(family.patch, this_family.patch) <= movement]
  my_helpers = [family
                for family in my_neighborhood
                if get_similarity(family.identity, this_family.identity) > cultural_distance ] # cultural_distance of my neighbors'
  # FIXME: YES, THE ORIGINAL SOURCE HAS A > THERE


def get_similarity(a, b):
  ap = a * weights_vector
  bp = b * weights_vector
  numerator = (ap * bp).sum()
  denominator = (ap * ap).sum() ** 0.5 * (bp * bp).sum() ** 0.5
  return numerator / denominator


def identify_groups():
  find_all_components
  for myself in families:
    if not myself.component > 0:
      continue
    my_group = [family for family in families if family.component == myself.component]
    group_size = len(my_group)


## The following two procedures are an adaptation from: Wilensky, U. (2005). NetLogo Giant Component model.
## http://ccl.northwestern.edu/netlogo/models/GiantComponent.
## Center for Connected Learning and Computer_Based Modeling, Northwestern University, Evanston, IL.
def find_all_components():
  for myself in families:
    myself.explored = false
    myself.group_size = 0
  for myself in families:
    if not link_neighbors():
      component = 0
      explored = True # families that don't cooperate (isolated agents) will have component = 0
  component_index = 0
  while True:
    if [family for family in families if not family.explored]:
        start = numpy.random.choice([family for family in families if not family.explored])
        component_index = ( component_index + 1)
        explore(start)
    else:
      break


## Finds all families reachable from this node
def explore(myself):
  if myself.explored:
    break
  myself.explored = true
  myself.component = component_index
  for neighbor in link_neighbors:
    explore(neighbor)


def decide_whether_to_move_or_stay():
  for myself in families:
    # 1. Calculate if I will be able to survive one more tick with my current level of energy (bearing in mind the depreciation factor) without working.
    # If that is the case, I leave with probability 0.05
    if (total_energy - survival_threshold) * conservation_factor > survival_threshold:
      if (numpy.random.random() < 0.05):
        move()
      # 2. I will have to work next tick because I will not be able to get by with my (depreciated) level of energy
      # Before moving, I will check if I will be able to get enough resources if I stay in the same patch
      resources_next_tick =  (myself.patch.max_resource - collected_energy) / 2 if (season == "hot") else myself.patch.max_resource
      if resources_next_tick * individual_capability  > survival_threshold:
        # If I can survive here, I will probably stay, but I will leave with probability 0.05
        if (numpy.random.random() < 0.05):
          move()
      else:
        move()


def move():
  if [patch for patch in patches if distance(patch, this_family.patch) <= movement if not patch.families]:
    number_of_movements = number_of_movements + 1
    move_to(numpy.random.choice([patch for patch in patches if distance(patch, this_family.patch) <= movement if not patch.families]))


def update_identity():
  # 1. Diffusion process (only for agents that have cooperated)
  if families:
    group_index = range(1, component_index + 1)
    for g in group_index:
      consensual_identity = compute_consensual_identity([family for family in families if family.component == g])
      for family in [family for family in families if family.component == g]
        if numpy.random.random() < 0.95:
          family.identity = consensual_identity
  # 2. Mutation process (for all families)
  for myself in families():
    if numpy.random.random() < internal_change_rate:
      index_vector = range(10)
      for i in index_vector:
          if numpy.random.random() >= (weights_vector[i]):
            myself.identity[i] = numpy.random.randint(1, 7)


def compute_consensual_identity(group):
  let consensus []
  foreach n_values 10 [?]
    [
      let consensual_trait one_of modes [ item ? identity] of group
      consensus = fput consensual_trait consensus
    ]
  consensus = reverse consensus
  report consensus


def update_technology
  if any? families
  [
    # 1. Diffusion process (only for agents that have cooperated)
    let group_index n_values component_index [1 + ?]
    foreach group_index
    [
      let this_group families with [component = ?]
      let average_technology_this_group mean [technology] of this_group
      ask this_group
      [
        ifelse technology < 0.95 * average_technology_this_group
        [
          if numpy.random.random() < 0.95
          [
            technology = technology + 0.1
            if technology > max [technology] of this_group [technology = max [technology] of this_group]
          ]
        ]
        [
          if technology < 1.05 * average_technology_this_group
            [
              if numpy.random.random() < 0.95
              [
                technology = technology + 0.01
                if technology > max [technology] of this_group [technology = max [technology] of this_group]
              ]
            ]
        ]
      ]
      ask this_group [conservation_factor = technology / 2]
    ]
  ]



  # 2. Mutation process (for all families)
  ask families
  [
    if numpy.random.random() < internal_change_rate
    [
      technology = random_normal average_technology diversity
      # possible correction to force technology to be in the interval [0 2] (and thus depreciation will be in [0 1]
      if technology < 0.02 [technology = 0.02]
      if technology > 2 [technology = 2]
      conservation_factor = technology / 2
    ]
  ]


def update_output_variables
  number_of_social_aggregates = component_index
  number_of_agents_in_social_aggregates = count families with [cooperation?]
  number_of_isolated_agents = count families with [not cooperation?]
  largest_group_size = ifelse_value any? (families with [group_size > 1]) [[group_size] of one_of families with_max [group_size]] ["N/A"]
  total_collected_energy = ifelse_value (any? families) [sum [collected_energy] of families] [0]
  collected_energy_standard_deviation = ifelse_value (count families >= 2) [standard_deviation [collected_energy] of families] [0]
  average_cultural_distance_in_aggregates = ifelse_value any? (families with [component > 0]) [mean [cultural_distance] of families with [component > 0]] ["N/A"]
  sd_cultural_distance_in_aggregates = ifelse_value any? (families with [component > 0]) [standard_deviation [cultural_distance] of families with [component > 0]] ["N/A"]
  total_number_of_starvation_deaths = total_number_of_starvation_deaths + number_of_agents_that_died_of_starvation
  sum_of_labor = sum [labor] of families
  if count families >= 2
  [
    mean_technology_of_families = mean [technology] of families
    std_technology_of_families = standard_deviation [technology] of families
  ]


def simulate():
  for step in range(1000):
    survive()
    update_resources_on_patches()
    update_cultural_distance()
    hunt_and_gather()
    identify_groups()
    update_identity()
    update_technology()
    decide_whether_to_move_or_stay()
    update_output_variables()
    tick()
    if not families:
      break


