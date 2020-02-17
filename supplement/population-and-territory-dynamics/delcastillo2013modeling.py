"""delcastillo2013modeling

Implement a geographic (abstracted) agent-based model of migrating
hunter-gatherer bands with culture. The implementation is based on the NetLogo
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

breed [families family]

globals
[
  season                         ; hot or cold
  weights-vector                 ; [0.5 0.25 0.12 0.06 0.03 0.015 0.0075 0.00375 0.001875 0.0009375]
  component-index                ; counter of groups

  #OUTPUT VARIABLES
  number-of-social-aggregates
  number-of-agents-in-social-aggregates
  number-of-isolated-agents
  number-of-agents-that-died-of-starvation
  number-of-living-agents      #(at the beginning of the tick)
  total-collected-energy
  collected-energy-standard-deviation
  average-cultural-distance-in-aggregates
  sd-cultural-distance-in-aggregates
  largest-group-size
  number-of-movements
  total-number-of-starvation-deaths
  sum-of-labor
  mean-technology-of-families
  std-technology-of-families
]


@attr.s
class Family:
  my_neighborhood = attr.ib()       # list of families around
  my_helpers = attr.ib()            # subset of my-neighborhood that will help me (with probability 95%)
  my_group = attr.ib()              # list of agents withins a social aggregate (i.e. families that belong to my group/component)
  group_size = attr.ib()            # size of the group (component) where I belong
  collected_energy = attr.ib()      # energy obtained during the current tick
  total_energy = attr.ib()          # total level of energy (including surplus)
  identity = attr.ib()              # 10-dimension vector with values 1 (no important) to 6 (very important)
  labor = attr.ib()                 # number of individuals within a household
  survival_threshold = attr.ib()    # 730 * labor
  technology = attr.ib()            # Gaussian-distributed: mean: average-technology; standard deviation: diversity
  conservation_factor = attr.ib()   # current value of technology divided by two
  individual_capability = attr.ib() # fraction of the resource that the agent is able to obtain when working on their own
  cultural_distance = attr.ib()     # cooperation will be posible if similarity >= cultural-distance
  cooperation = attr.ib()           # true if I have been involved in a cooperation process (i.e. I have helped someone or someone has helped me)
  explored = attr.ib()              # needed for detecting components
  component = attr.ib()             # group (i.e. component in the network) where the family belongs


@attr.s
class Patch:
  difficulty = attr.ib()      # (h)
  max_resource = attr.ib()    # maximum amount of resource (in hot season)
  resource = attr.ib()        # current amount of resource
  amount-consumed = attr.ib() # amount of energy consumed during the current tick
  exploited = attr.ib()       # true during one tick after been exploited


def setup():
  clear_all()
  reset_ticks()

  #SET GLOBAL VARIABLES
  weights_vector = [0.5 0.25 0.12 0.06 0.03 0.015 0.0075 0.00375 0.001875 0.0009375]

  #CREATE AND INITIALIZE FAMILIES
  initial_identity = numpy.random.randint(1, 7)  #all the families start with the same identity vector

  ask n-of initial-population patches
  [
    sprout-families 1
    [
      set shape "person"
      set color red
      set identity initial-identity
      set labor random-poisson labor-average if labor < 2 [set labor 2] #(correction for cases with labor < 2)
      set survival-threshold (730 * labor) / 2
      set technology random-normal average-technology diversity
          # correction to force technology to be in the interval [0 2] (and thus depreciation will be in [0 1]
          if technology < 0.02 [set technology 0.02]
          if technology > 2 [set technology 2]
      set conservation-factor technology / 2
      set total-energy 1 + survival-threshold * (1 + conservation-factor) / conservation-factor    # this means that all the families will have enough energy to survive the first tick without hunting
      set collected-energy 0
      set my-group families with [who < 0]     # initialize my-group = nobody
    ]
  ]


  #INITIALIZE PATCHES
  ask patches
  [
    # Resources are uniformly distributed in [max-resource-on-patches / 1000, max-resource-on-patches]
    set max-resource (max-resource-on-patches / 1000) + random ( max-resource-on-patches - (max-resource-on-patches / 1000) )
    set difficulty random-float 1
    set exploited? false
    set amount-consumed 0
    set pcolor 69
  ]


def go():
  survive
  update-resources-on-patches
  update-cultural-distance
  hunt-and-gather
  identify-groups
  update-identity
  update-technology
  decide-whether-to-move-or-stay
  update-output-variables
  tick
  #if not any? families [stop]
  if ticks = 1000 or not any? families [stop]


def survive():
  set number-of-agents-that-died-of-starvation 0
  set number-of-living-agents count families

  # 1. Substract survival-threshold. If they didn't reach their survival threshold in the previous tick, their current level of total-energy will be < 0 at this point;
  # in this case their amount of labor is reduced by one unit
  ask families
    [
      set total-energy total-energy - survival-threshold

      if total-energy < 0
      [
        set labor labor - 1
        set survival-threshold (730 * labor) / 2
        set total-energy 0

        #if their amount of labor goes below 2, the family (agent) dies
        if labor < 2
        [
          set number-of-agents-that-died-of-starvation number-of-agents-that-died-of-starvation + 1
          die
        ]
      ]
    ]

  # 2. Vegatative reproduction: The amount of labor within a family goes up one unit every 6 ticks
  ask families
  [
    if ticks > 1 and remainder ticks 30 = 0
    [
      set labor labor + 1
      set survival-threshold (730 * labor) / 2
    ]
  ]

  # if the amount of labor reaches the value 10, the family is split into two families with probability 95%
  # (provided there is enough room in the neighbourhood for another family: there cannot be more than one family in the same patch)
  ask families
  [
    if labor >= 10
    [
      if random-float 1 < .95
      [
        if any? patches in-radius movement with [not any? turtles-here]
          [
            set labor floor (labor / 2)
            set survival-threshold (730 * labor) / 2
            hatch-families 1 [move-to one-of patches in-radius movement with [not any? turtles-here]]
          ]
      ]
    ]
  ]


def update-resources-on-patches():
  set season ifelse-value (remainder ticks  2 = 0) ["hot"] ["cold"]
  ask patches
  [
    ifelse exploited?
    [
      # patches that were exploited in the previous tick:
      # if the current season is cold, the patch will not regenerate completely: (max-resource - amount-consumed)/2
      # if the current season is hot, the level of resource will be max-resource (no matter whether it was exploited during the previous tick or not)
      if season = "hot" [set resource max-resource]
      if season = "cold" [set resource (max-resource - amount-consumed) / 2]
    ]
    [
      set resource ifelse-value (season = "hot") [max-resource] [max-resource / 2]
    ]
    set exploited? false
    set amount-consumed 0
  ]
  ifelse (display-resorces?) [ask patches [set pcolor scale-color green resource 100000 100]] [ask patches [set pcolor 69]]


def update-cultural-distance:
  ask families
  [
    ifelse ([resource] of patch-here) <= survival-threshold
    [
      # They will refuse to cooperate, as they are unable to reach their survival-threshold no matter the number of people who help them
      set cultural-distance 1
    ]
    [
      let extra-workers max list 0 (  (survival-threshold / (difficulty * (([resource] of patch-here) - survival-threshold) ) ) ^ ( 1 / technology) - labor)
      # extra-workers could be negative (if I have plenty of labor to obtain my survival threshold: I don't need any help from anyone)  --> in this case extra-workers = 0
      set cultural-distance min list 1 ((1 / 100) * extra-workers) # cultural distance could be greater than 1 if the number of extra workers goes above 100 --> --> in this case cultural distance = 1
    ]
  ]


def hunt-and-gather():
  ask families
  [
    set total-energy total-energy * conservation-factor #depreciate agents' energy at the beginning of the tick
    set cooperation? false
    set my-helpers families with [who < 0] # (empty agentlist)
    set my-group families with [who < 0]   # (empty agentlist)
    set collected-energy 0
  ]

  ask links [die]


  # families with total-energy > survival-threshold don't hunt or gather

  ask families with [total-energy < survival-threshold]
  [
    # Try to act individually:
    set individual-capability calculate-individual-capability
    let productivity ([resource] of patch-here) * individual-capability

    ifelse cooperation-allowed?


    [
      # FAMILIES ARE ALLOWED TO ASK FOR HELP WHENEVER THEY ARE UNABLE TO REACH THEIR SURVIVAL THRESHOLD
      ifelse productivity >= survival-threshold
      [
        # ACT INDIVIDUALLY
        # I don't need to cooperate: I will act individually and collect as much energy as I need to reach my survival-threshold

        set collected-energy min list (survival-threshold - total-energy) (individual-capability) * [resource] of patch-here
        set total-energy total-energy + collected-energy
        ask patch-here #update resource on patch:
          [
            set amount-consumed [collected-energy] of myself
            set resource resource - amount-consumed
            set exploited? true
          ]
      ]
      [
        # COOPERATION (I need to ask for help)
        let capability individual-capability
        identify-neighbors # define my-helpers (my-helpers contains a list of families that can potentially help me, i.e. our similarity is greater than their cultural-distance)
        let agents-willing-to-help ask-for-help
        if length agents-willing-to-help > 0 # if someone helps me, my capability will be aggregated-capability. Otherwise it will be my individual-capability
        [
          set capability calculate-aggregated-capability agents-willing-to-help
        ]

        set collected-energy ([resource] of patch-here) * capability

        set total-energy total-energy + collected-energy        #(therefore, total-energy might be greater than her survival-threshold

        #update resource on patch:
        ask patch-here
          [
            set amount-consumed [collected-energy] of myself
            set resource resource - amount-consumed
            set exploited? true
          ]
      ]
    ]



    [
      # FAMILIES ARE NOT ALLOWED TO ASK FOR HELP IF THEY ARE UNABLE TO REACH THEIR SURVIVAL THRESHOLD. THEY WILL ACT INDIVIDUALLY

      set collected-energy min list (survival-threshold - total-energy) (individual-capability) * [resource] of patch-here
      set total-energy total-energy + collected-energy
      ask patch-here #update resource on patch:
          [
            set amount-consumed [collected-energy] of myself
            set resource resource - amount-consumed
            set exploited? true
          ]
    ]
  ]


def -report calculate-individual-capability
  report 1 / (1 + (1 / (([difficulty] of patch-here) * (labor ^ technology))))


def -report ask-for-help
  # Determine the list of agents that are willing to cooperate
  let agents-willing-to-cooperate []
  if count my-helpers > 0  # if my list of helpers is not empty
  [
    #each one of my helpers will help me with p=95%
    foreach [self] of my-helpers
    [
      if random-float 1 < 0.95
      [
        set agents-willing-to-cooperate fput ? agents-willing-to-cooperate
      ]
    ]
    # The agents contained in the list agents-willing-to-cooperate will help me
    ask turtle-set agents-willing-to-cooperate [create-link-with myself]   #helpers create a link with helped agents
    ask links [set color black set thickness .1]
    set cooperation? true # I have cooperated...
    ask turtle-set agents-willing-to-cooperate [set cooperation? true]  # ... and so have my helpers
  ]
  report agents-willing-to-cooperate


def -report calculate-aggregated-capability [agents-willing-to-cooperate]
  let max-technology-within-my-local-group max [technology] of (turtle-set agents-willing-to-cooperate self)   # maximum technology within my group (including myself)
  let returns-to-cooperation 2 - ([resource] of patch-here / max-resource-on-patches)
  let total-labor labor + sum [labor] of turtle-set agents-willing-to-cooperate                                  # total labor within my group (including myself)
  let aggregated-capability 1 / (1 + ( 1 / ( ( ([difficulty] of patch-here) * (total-labor ^ max-technology-within-my-local-group)) ^ returns-to-cooperation ) ) )
  report aggregated-capability


def identify-neighbors
  set my-neighborhood other families-on patches in-radius movement
  set my-helpers my-neighborhood with [(get-similarity identity ([identity] of myself)) > cultural-distance ] # cultural-distance of my neighbors'


def -report get-similarity [a b]
  let ap (map [?1 * ?2] a weights-vector)
  let bp (map [?1 * ?2] b weights-vector)
  let numerator sum (map [?1 * ?2] ap bp)
  let denominator sqrt((sum map [? * ?] ap)) * sqrt((sum map [? * ?] bp))
  report numerator / denominator


def identify-groups
  find-all-components
  ask families with [component > 0]
  [
    set my-group families with [component = [component] of myself]
    set group-size count my-group
  ]


## The following two procedures are an adaptation from: Wilensky, U. (2005). NetLogo Giant Component model.
## http://ccl.northwestern.edu/netlogo/models/GiantComponent.
## Center for Connected Learning and Computer-Based Modeling, Northwestern University, Evanston, IL.
def find-all-components
  ask families [set explored? false set group-size 0]
  ask families with [not any? link-neighbors] [set component 0 set explored? true] # families that don't cooperate (isolated agents) will have component = 0
  set component-index 0
  loop
  [ ifelse any? families with [not explored?]
      [
        let start one-of sort families with [ not explored? ]
        set component-index ( component-index + 1)
        ask start [ explore ]
      ]
      [stop]
  ]


## Finds all families reachable from this node
def explore
  if explored? [stop]
  set explored? true
  set component component-index
  ask link-neighbors [explore]


def decide-whether-to-move-or-stay
  ask families
  [
    # 1. Calculate if I will be able to survive one more tick with my current level of energy (bearing in mind the depreciation factor) without working.
    # If that is the case, I leave with probability 0.05
    ifelse (total-energy - survival-threshold) * conservation-factor > survival-threshold
    [
      if (random-float 1 < 0.05) [move]
    ]
    [
      # 2. I will have to work next tick because I will not be able to get by with my (depreciated) level of energy
      # Before moving, I will check if I will be able to get enough resources if I stay in the same patch
      let resources-next-tick ifelse-value (season = "hot") [([max-resource] of patch-here - collected-energy) / 2] [[max-resource] of patch-here]
      ifelse resources-next-tick * individual-capability  > survival-threshold
      [
        # If I can survive here, I will probably stay, but I will leave with probability 0.05
        if (random-float 1 < 0.05) [move]
      ]
      [
        move
      ]
    ]
  ]


def move
  if any? patches in-radius movement with [not any? turtles-here]
  [
    set number-of-movements number-of-movements + 1
    move-to one-of patches in-radius movement with [not any? turtles-here]
  ]


def update-identity
  # 1. Diffusion process (only for agents that have cooperated)
  if any? families
  [
    let group-index n-values component-index [1 + ?]
    foreach group-index
    [
      let consensual-identity compute-consensual-identity families with [component = ?]
      ask families with [component = ?]
      [
        if random-float 1 < 0.95 [set identity consensual-identity]
      ]
    ]
  ]
  # 2. Mutation process (for all families)
  ask families
  [
    if random-float 1 < internal-change-rate
    [
      let index-vector n-values 10 [?]
      foreach index-vector
        [
          if random-float 1 >= (item ? weights-vector) [set identity replace-item ? identity (1 + random 5)]
        ]
    ]
  ]


def -report compute-consensual-identity [group]
  let consensus []
  foreach n-values 10 [?]
    [
      let consensual-trait one-of modes [ item ? identity] of group
      set consensus fput consensual-trait consensus
    ]
  set consensus reverse consensus
  report consensus


def update-technology
  if any? families
  [
    # 1. Diffusion process (only for agents that have cooperated)
    let group-index n-values component-index [1 + ?]
    foreach group-index
    [
      let this-group families with [component = ?]
      let average-technology-this-group mean [technology] of this-group
      ask this-group
      [
        ifelse technology < 0.95 * average-technology-this-group
        [
          if random-float 1 < 0.95
          [
            set technology technology + 0.1
            if technology > max [technology] of this-group [set technology max [technology] of this-group]
          ]
        ]
        [
          if technology < 1.05 * average-technology-this-group
            [
              if random-float 1 < 0.95
              [
                set technology technology + 0.01
                if technology > max [technology] of this-group [set technology max [technology] of this-group]
              ]
            ]
        ]
      ]
      ask this-group [set conservation-factor technology / 2]
    ]
  ]



  # 2. Mutation process (for all families)
  ask families
  [
    if random-float 1 < internal-change-rate
    [
      set technology random-normal average-technology diversity
      # possible correction to force technology to be in the interval [0 2] (and thus depreciation will be in [0 1]
      if technology < 0.02 [set technology 0.02]
      if technology > 2 [set technology 2]
      set conservation-factor technology / 2
    ]
  ]


def update-output-variables
  set number-of-social-aggregates component-index
  set number-of-agents-in-social-aggregates count families with [cooperation?]
  set number-of-isolated-agents count families with [not cooperation?]
  set largest-group-size ifelse-value any? (families with [group-size > 1]) [[group-size] of one-of families with-max [group-size]] ["N/A"]
  set total-collected-energy ifelse-value (any? families) [sum [collected-energy] of families] [0]
  set collected-energy-standard-deviation ifelse-value (count families >= 2) [standard-deviation [collected-energy] of families] [0]
  set average-cultural-distance-in-aggregates ifelse-value any? (families with [component > 0]) [mean [cultural-distance] of families with [component > 0]] ["N/A"]
  set sd-cultural-distance-in-aggregates ifelse-value any? (families with [component > 0]) [standard-deviation [cultural-distance] of families with [component > 0]] ["N/A"]
  set total-number-of-starvation-deaths total-number-of-starvation-deaths + number-of-agents-that-died-of-starvation
  set sum-of-labor sum [labor] of families
  if count families >= 2
  [
    set mean-technology-of-families mean [technology] of families
    set std-technology-of-families  standard-deviation [technology] of families
  ]
