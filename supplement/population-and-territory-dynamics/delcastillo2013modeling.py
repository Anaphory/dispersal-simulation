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
  weights_vector = [0.5 0.25 0.12 0.06 0.03 0.015 0.0075 0.00375 0.001875 0.0009375]

  #CREATE AND INITIALIZE FAMILIES
  initial_identity = numpy.random.randint(1, 7)  #all the families start with the same identity vector

  ask n_of initial_population patches
  [
    sprout_families 1
    [
      shape = "person"
      color = red
      identity = initial_identity
      labor = random_poisson labor_average if labor < 2 [set labor 2] #(correction for cases with labor < 2)
      survival_threshold = (730 * labor) / 2
      technology = random_normal average_technology diversity
          # correction to force technology to be in the interval [0 2] (and thus depreciation will be in [0 1]
          if technology < 0.02 [technology = 0.02]
          if technology > 2 [technology = 2]
      conservation_factor = technology / 2
      total_energy = 1 + survival_threshold * (1 + conservation_factor) / conservation_factor    # this means that all the families will have enough energy to survive the first tick without hunting
      collected_energy = 0
      my_group = families with [who < 0]     # initialize my_group = nobody
    ]
  ]


  #INITIALIZE PATCHES
  ask patches
  [
    # Resources are uniformly distributed in [max_resource_on_patches / 1000, max_resource_on_patches]
    max_resource = (max_resource_on_patches / 1000) + random ( max_resource_on_patches - (max_resource_on_patches / 1000) )
    difficulty = random_float 1
    exploited? = false
    amount_consumed = 0
    pcolor = 69
  ]


def go():
  survive
  update_resources_on_patches
  update_cultural_distance
  hunt_and_gather
  identify_groups
  update_identity
  update_technology
  decide_whether_to_move_or_stay
  update_output_variables
  tick
  #if not any? families [stop]
  if ticks = 1000 or not any? families [stop]


def survive():
  number_of_agents_that_died_of_starvation = 0
  number_of_living_agents = count families

  # 1. Substract survival_threshold. If they didn't reach their survival threshold in the previous tick, their current level of total_energy will be < 0 at this point;
  # in this case their amount of labor is reduced by one unit
  ask families
    [
      total_energy = total_energy - survival_threshold

      if total_energy < 0
      [
        labor = labor - 1
        survival_threshold = (730 * labor) / 2
        total_energy = 0

        #if their amount of labor goes below 2, the family (agent) dies
        if labor < 2
        [
          number_of_agents_that_died_of_starvation = number_of_agents_that_died_of_starvation + 1
          die
        ]
      ]
    ]

  # 2. Vegatative reproduction: The amount of labor within a family goes up one unit every 6 ticks
  ask families
  [
    if ticks > 1 and remainder ticks 30 = 0
    [
      labor = labor + 1
      survival_threshold = (730 * labor) / 2
    ]
  ]

  # if the amount of labor reaches the value 10, the family is split into two families with probability 95%
  # (provided there is enough room in the neighbourhood for another family: there cannot be more than one family in the same patch)
  ask families
  [
    if labor >= 10
    [
      if random_float 1 < .95
      [
        if any? patches in_radius movement with [not any? turtles_here]
          [
            labor = floor (labor / 2)
            survival_threshold = (730 * labor) / 2
            hatch_families 1 [move_to one_of patches in_radius movement with [not any? turtles_here]]
          ]
      ]
    ]
  ]


def update_resources_on_patches():
  season = ifelse_value (remainder ticks  2 = 0) ["hot"] ["cold"]
  ask patches
  [
    ifelse exploited?
    [
      # patches that were exploited in the previous tick:
      # if the current season is cold, the patch will not regenerate completely: (max_resource - amount_consumed)/2
      # if the current season is hot, the level of resource will be max_resource (no matter whether it was exploited during the previous tick or not)
      if season = "hot" [resource = max_resource]
      if season = "cold" [resource = (max_resource - amount_consumed) / 2]
    ]
    [
      resource = ifelse_value (season = "hot") [max_resource] [max_resource / 2]
    ]
    exploited? = false
    amount_consumed = 0
  ]
  ifelse (display_resorces?) [ask patches [pcolor = scale_color green resource 100000 100]] [ask patches [set pcolor 69]]


def update_cultural_distance:
  ask families
  [
    ifelse ([resource] of patch_here) <= survival_threshold
    [
      # They will refuse to cooperate, as they are unable to reach their survival_threshold no matter the number of people who help them
      cultural_distance = 1
    ]
    [
      let extra_workers max list 0 (  (survival_threshold / (difficulty * (([resource] of patch_here) - survival_threshold) ) ) ^ ( 1 / technology) - labor)
      # extra_workers could be negative (if I have plenty of labor to obtain my survival threshold: I don't need any help from anyone)  --> in this case extra_workers = 0
      cultural_distance = min list 1 ((1 / 100) * extra_workers) # cultural distance could be greater than 1 if the number of extra workers goes above 100 --> --> in this case cultural distance = 1
    ]
  ]


def hunt_and_gather():
  ask families
  [
    total_energy = total_energy * conservation_factor #depreciate agents' energy at the beginning of the tick
    cooperation? = false
    my_helpers = families with [who < 0] # (empty agentlist)
    my_group = families with [who < 0]   # (empty agentlist)
    collected_energy = 0
  ]

  ask links [die]


  # families with total_energy > survival_threshold don't hunt or gather

  ask families with [total_energy < survival_threshold]
  [
    # Try to act individually:
    individual_capability = calculate_individual_capability
    let productivity ([resource] of patch_here) * individual_capability

    ifelse cooperation_allowed?


    [
      # FAMILIES ARE ALLOWED TO ASK FOR HELP WHENEVER THEY ARE UNABLE TO REACH THEIR SURVIVAL THRESHOLD
      ifelse productivity >= survival_threshold
      [
        # ACT INDIVIDUALLY
        # I don't need to cooperate: I will act individually and collect as much energy as I need to reach my survival_threshold

        collected_energy = min list (survival_threshold - total_energy) (individual_capability) * [resource] of patch_here
        total_energy = total_energy + collected_energy
        ask patch_here #update resource on patch:
          [
            amount_consumed = [collected_energy] of myself
            resource = resource - amount_consumed
            exploited? = true
          ]
      ]
      [
        # COOPERATION (I need to ask for help)
        let capability individual_capability
        identify_neighbors # define my_helpers (my_helpers contains a list of families that can potentially help me, i.e. our similarity is greater than their cultural_distance)
        let agents_willing_to_help ask_for_help
        if length agents_willing_to_help > 0 # if someone helps me, my capability will be aggregated_capability. Otherwise it will be my individual_capability
        [
          capability = calculate_aggregated_capability agents_willing_to_help
        ]

        collected_energy = ([resource] of patch_here) * capability

        total_energy = total_energy + collected_energy        #(therefore, total_energy might be greater than her survival_threshold

        #update resource on patch:
        ask patch_here
          [
            amount_consumed = [collected_energy] of myself
            resource = resource - amount_consumed
            exploited? = true
          ]
      ]
    ]



    [
      # FAMILIES ARE NOT ALLOWED TO ASK FOR HELP IF THEY ARE UNABLE TO REACH THEIR SURVIVAL THRESHOLD. THEY WILL ACT INDIVIDUALLY

      collected_energy = min list (survival_threshold - total_energy) (individual_capability) * [resource] of patch_here
      total_energy = total_energy + collected_energy
      ask patch_here #update resource on patch:
          [
            amount_consumed = [collected_energy] of myself
            resource = resource - amount_consumed
            exploited? = true
          ]
    ]
  ]


def -report calculate_individual_capability
  report 1 / (1 + (1 / (([difficulty] of patch_here) * (labor ^ technology))))


def -report ask_for_help
  # Determine the list of agents that are willing to cooperate
  let agents_willing_to_cooperate []
  if count my_helpers > 0  # if my list of helpers is not empty
  [
    #each one of my helpers will help me with p=95%
    foreach [self] of my_helpers
    [
      if random_float 1 < 0.95
      [
        agents_willing_to_cooperate = fput ? agents_willing_to_cooperate
      ]
    ]
    # The agents contained in the list agents_willing_to_cooperate will help me
    ask turtle_agents_willing_to_cooperate = [create_link_with myself]   #helpers create a link with helped agents
    ask links [color = black set thickness .1]
    cooperation? = true # I have cooperated...
    ask turtle_agents_willing_to_cooperate = [set cooperation? true]  # ... and so have my helpers
  ]
  report agents_willing_to_cooperate


def -report calculate_aggregated_capability [agents_willing_to_cooperate]
  let max_technology_within_my_local_group max [technology] of (turtle_agents_willing_to_cooperate = self)   # maximum technology within my group (including myself)
  let returns_to_cooperation 2 - ([resource] of patch_here / max_resource_on_patches)
  let total_labor labor + sum [labor] of turtle_agents_willing_to_cooperate = # total labor within my group (including myself)
  let aggregated_capability 1 / (1 + ( 1 / ( ( ([difficulty] of patch_here) * (total_labor ^ max_technology_within_my_local_group)) ^ returns_to_cooperation ) ) )
  report aggregated_capability


def identify_neighbors
  my_neighborhood = other families_on patches in_radius movement
  my_helpers = my_neighborhood with [(get_similarity identity ([identity] of myself)) > cultural_distance ] # cultural_distance of my neighbors'


def -report get_similarity [a b]
  let ap (map [?1 * ?2] a weights_vector)
  let bp (map [?1 * ?2] b weights_vector)
  let numerator sum (map [?1 * ?2] ap bp)
  let denominator sqrt((sum map [? * ?] ap)) * sqrt((sum map [? * ?] bp))
  report numerator / denominator


def identify_groups
  find_all_components
  ask families with [component > 0]
  [
    my_group = families with [component = [component] of myself]
    group_size = count my_group
  ]


## The following two procedures are an adaptation from: Wilensky, U. (2005). NetLogo Giant Component model.
## http://ccl.northwestern.edu/netlogo/models/GiantComponent.
## Center for Connected Learning and Computer_Based Modeling, Northwestern University, Evanston, IL.
def find_all_components
  ask families [explored? = false set group_size 0]
  ask families with [not any? link_neighbors] [component = 0 set explored? true] # families that don't cooperate (isolated agents) will have component = 0
  component_index = 0
  loop
  [ ifelse any? families with [not explored?]
      [
        let start one_of sort families with [ not explored? ]
        component_index = ( component_index + 1)
        ask start [ explore ]
      ]
      [stop]
  ]


## Finds all families reachable from this node
def explore
  if explored? [stop]
  explored? = true
  component = component_index
  ask link_neighbors [explore]


def decide_whether_to_move_or_stay
  ask families
  [
    # 1. Calculate if I will be able to survive one more tick with my current level of energy (bearing in mind the depreciation factor) without working.
    # If that is the case, I leave with probability 0.05
    ifelse (total_energy - survival_threshold) * conservation_factor > survival_threshold
    [
      if (random_float 1 < 0.05) [move]
    ]
    [
      # 2. I will have to work next tick because I will not be able to get by with my (depreciated) level of energy
      # Before moving, I will check if I will be able to get enough resources if I stay in the same patch
      let resources_next_tick ifelse_value (season = "hot") [([max_resource] of patch_here - collected_energy) / 2] [[max_resource] of patch_here]
      ifelse resources_next_tick * individual_capability  > survival_threshold
      [
        # If I can survive here, I will probably stay, but I will leave with probability 0.05
        if (random_float 1 < 0.05) [move]
      ]
      [
        move
      ]
    ]
  ]


def move
  if any? patches in_radius movement with [not any? turtles_here]
  [
    number_of_movements = number_of_movements + 1
    move_to one_of patches in_radius movement with [not any? turtles_here]
  ]


def update_identity
  # 1. Diffusion process (only for agents that have cooperated)
  if any? families
  [
    let group_index n_values component_index [1 + ?]
    foreach group_index
    [
      let consensual_identity compute_consensual_identity families with [component = ?]
      ask families with [component = ?]
      [
        if random_float 1 < 0.95 [identity = consensual_identity]
      ]
    ]
  ]
  # 2. Mutation process (for all families)
  ask families
  [
    if random_float 1 < internal_change_rate
    [
      let index_vector n_values 10 [?]
      foreach index_vector
        [
          if random_float 1 >= (item ? weights_vector) [identity = replace_item ? identity (1 + random 5)]
        ]
    ]
  ]


def -report compute_consensual_identity [group]
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
          if random_float 1 < 0.95
          [
            technology = technology + 0.1
            if technology > max [technology] of this_group [technology = max [technology] of this_group]
          ]
        ]
        [
          if technology < 1.05 * average_technology_this_group
            [
              if random_float 1 < 0.95
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
    if random_float 1 < internal_change_rate
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
