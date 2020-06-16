use petgraph;
use petgraph::algo::dijkstra;

// https://www.sciencedirect.com/science/article/pii/S0305440307000635?casa_token=9pFpl8kYoGIAAAAA:jjKyuCgnOATkREq3jeZjI9mTlbZMZLQHQUi7mL-kej7Z97xA6EBzEcUWHXPUMWux_sgkdkesOw
// Specifically, the time invested in caching acorn in dispersed locations should comprise no more than a day's labor to build, fill, and later re-access cached acorn, this again based on the assumption that hunter–gatherers tend to work no more than 8 h a day. In its simplest sense then, caches should be no more than a half-day round trip from winter settlements, with the other half-day left for other activities. Using this logic, a central place cacher walking 4.7 km/h travels 18.8 km round trip in half a workday (4 h), or 9.4 km one-way, this latter figure the predicted foraging limit encircling winter settlements.

// iwamura2014agentbaseda: Walking speed in grassland 4km/h; Walking speed in forest 2km/h; BMR 279.6 kJ/h

// For each movement from one cell to another, a household agent loses energy, c_walk:
// $c_walk = 3.2 β h_walk$
// β represents the energy consumption of 1 BMR in Joules, h_walk is the number of hours required to walk through a cell. h walk depends on the size of cells and the walking speed, assumed to be 3 km per hour.

/**
Base metabolic rate in kcal/h

Following [@harris1918biometric], the base metabolic rate (in kcal/day) can be
estimated from the age in *years*, the height in *cm*, and the weight in *kg*,
using different linear equations for men and women.

*/
pub fn bmr(age: f32, height: f32, weight: f32, male: bool) -> f32 {
    (if male {
        66. + (13.7 * weight) + (5. * height) - (6.8 * age)
    } else {
        655. + (9.6 * weight) + (1.8 * height) - (4.7 * age)
    }) / 24. // kcal/day → kcal/h
}

/**
We assume the base metabolic rate for an average 30-year old male of 170 cm
height and 65 kg weight. For this configuration, [@iwamura2014agentbased] cite a
BMR of 279.6 kJ/h, which matches our `bmr` function up to a rounding error.

```
extern crate approx;
use approx::*;

# use minimal_metabolic_cost_surface::{BMR, bmr};

assert_abs_diff_eq!(
    BMR,
    bmr(30., 170., 65., true));

assert_abs_diff_eq!(
    279.6,

    BMR *
    4.184,  // kJ / kcal
epsilon=0.5);
```

 */
pub static BMR: f32 = 66.770836; // kcal/h

/**
Compute the metabolic rate for walking, in kcal/h

[@wood2006energetically] estimate the metabolic rate for an individual walking
through non-flat terrain. They provide the results in Watts. For convenience of
compatibility with our other functions, which are computing energy in kcal,
distances in km and velocities in km/h, we compute the result instead in kcal/h.
*/
pub fn metabolic_rate(
    // in kg
    subject_weight: f32,
    // in kg
    load: f32,
    // in km/h
    walking_speed: f32,
    // in % (m/100m)
    grade: f32,
    // 1 for treadmill walking
    terrain_factor: f32,
    // standing metabolic rate, the minimal metabolic rate which walking should never be less than
    standing_metabolic_rate: f32
) -> f32 {
    let w = subject_weight;
    let l = load;
    let v = walking_speed / 3.6; // now in m/s as required by the reference formula
    let g = grade;
    let eta = terrain_factor;

    let m = 1.5 * w + 2.0 * (w + l) / (l / w).powi(2) + eta * (w + l) * (1.5 * v.powi(2) + 0.35 * v * g);
    let c = if g < 0. {
        eta * ((g * (w + l) * v) / 3.5 - (((w + l) * (g + 6.).powi(2)) / w) + (25. - v.powi(2)))
    } else {
        0.0
    };

    f32::max((m - c) * (3600. / 4184.), // W → kcal/h
             standing_metabolic_rate)
}

pub fn sqrt(w: f32) -> f32 {
    w.powf(0.5)
}

/**
Compute the energetically optimal speed for some gradient

```
use approx::*;

# use minimal_metabolic_cost_surface::{optimal_speed,metabolic_rate,SMR};

let v = optimal_speed(65., 10., 0.0001, 1., SMR);

assert!(metabolic_rate(65., 10., v, 0.0001, 1., SMR)/v <  metabolic_rate(65., 10., v+1., 0.0001, 1., SMR) / (v+1.));
assert!(metabolic_rate(65., 10., v, 0.0001, 1., SMR)/v <  metabolic_rate(65., 10., v-1., 0.0001, 1., SMR) / (v-1.));

assert_abs_diff_eq!(
5.0,
v);
```

*/
pub fn optimal_speed(
    // in kg
    subject_weight: f32,
    // in kg
    load: f32,
    // in % (m/100m)
    grade: f32,
    // 1 for treadmill walking
    terrain_factor: f32,
    // standing metabolic rate, the minimal metabolic rate which walking should never be less than
    standing_metabolic_rate: f32
) -> f32 {
    let w = subject_weight;
    let l = load;
    let g = grade;
    let eta = terrain_factor;
    0.0 * 3.6
}

pub static SMR: f32 = 1.2 * BMR;
pub static SPEED: f32 = 4.; //km/h

/// in kcal
pub fn travel_cost(
    // in km
    horizontal_distance: f32,
    // destination elevation minus starting elevation, in m
    ascent: f32,
) -> f32 {
    metabolic_rate(65., 10., SPEED, ascent / horizontal_distance * 0.1, 1.0, SMR) * horizontal_distance / SPEED
}

/**
Using the slope in %, calculate the navigation speed in m/s

This function calculates the off-road navigation speed (for male cadets in
forested areas with navigational targets every 20 minutes) following
[@irmischer2018measuring].
 */
pub fn travel_speed(
    // in %
    slope: f32
) -> f32 {
    return 0.11 + 0.67 * f32::exp(-(slope + 2.0).powi(2) / 1800.)
}
