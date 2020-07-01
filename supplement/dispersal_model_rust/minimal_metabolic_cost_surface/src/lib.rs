use petgraph;
use std::f64::consts::PI;


pub mod hexgrid;

pub type KCal = f32;

// TODO: Maybe KCal should be a finite f32? If resources become infinite,
// something has gone wrong. I bet there is a Rust crate for that…

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

pub static SMR: f32 = 1.2 * BMR;

/// in kcal
pub fn navigation_cost(
    // in km
    horizontal_distance: f32,
    // destination elevation minus starting elevation, in m
    ascent: f32,
) -> f32 {
    let slope = ascent / horizontal_distance * 0.1; // in %
    let speed = navigation_speed(slope) * 3.6; // m/s → km/h
    metabolic_rate(65., 10., speed, slope, 1.0, SMR) * horizontal_distance / speed
}

/**
Using the slope in %, calculate the navigation speed in m/s

This function calculates the off-road navigation speed (for male cadets in
forested areas with navigational targets every 20 minutes) following
[@irmischer2018measuring].

> [T]he fastest off-road navigation speed was 0.78 m/s for males […] with a peak
> at −2%.
``` use approx::*;

# use minimal_metabolic_cost_surface::navigation_speed;

assert_abs_diff_eq!(
  navigation_speed(-2.),
  0.78);

assert!(navigation_speed(-2.) > navigation_speed(-1.5))
assert!(navigation_speed(-2.) > navigation_speed(-2.5))
```

 */
pub fn navigation_speed(
    // in %
    slope: f32
) -> f32 {
    return 0.11 + 0.67 * f32::exp(-(slope + 2.0).powi(2) / 1800.)
}

pub struct Coords {
    x: f32,
    y: f32
}

pub fn navigation_distance(
    start_coordinates: Coords,
    end_coordinates: Vec<Coords>,
    dem: f32
) -> Vec<f32> {
    vec!()
}


pub fn indices_from_coordinates(
    coordinates: hexgrid::GeoCoord,
    width: usize,
    total: usize
) -> (usize, usize) {
    let column = ((coordinates.lon + PI) / PI / 2. * width as f64).round() as usize;
    let row = ((-coordinates.lat + PI / 2.) / PI * (total / width) as f64)
        .round() as usize;
    (row, column)
}


pub fn data_from_coordinates<'a>(
    coordinates: hexgrid::GeoCoord,
    image_pixels: &'a Vec<u32>,
    pixels_width: usize,
) -> Option<&'a u32> {
    let (row, column) = indices_from_coordinates(coordinates, pixels_width, image_pixels.len());
    let index = row * pixels_width + column;
    image_pixels.get(index)
}


pub fn load_tif(path: &std::path::Path) -> Option<(Vec<u32>, usize)> {
    // TODO: How do I specify data paths?
    let f = std::fs::File::open(path).ok()?;
    let maybe_image = tiff::decoder::Decoder::new(f);
    let mut image = maybe_image.ok()?;
    let (width, height) = image.dimensions().ok()?;
    println!("# {}x{}", width, height);
    let outcome = image.read_image().ok()?;
    println!("# Image read");
    let vec = match outcome {
        tiff::decoder::DecodingResult::U8(v) => v.iter().map(|g| u32::from(*g)).collect(),
        tiff::decoder::DecodingResult::U16(v) => v.iter().map(|g| u32::from(*g)).collect(),
        tiff::decoder::DecodingResult::U32(w) => w,
        tiff::decoder::DecodingResult::U64(v) => v.iter().map(|g| *g as u32).collect(),
    };
    Some((vec, width as usize))
}

use std::path::Path;

/**
Compute the pairwise, directed costs to move from one hex to its neighbors
*/
pub fn pairwise_costs(candidates: Vec<hexgrid::Index>) -> Option<HashMap<hexgrid::Index, HashMap<hexgrid::Index, KCal>>> {
    let (pixels, width) = load_tif(Path::new("dem.tif"))?;
    let total = pixels.len();
    candidates
        .iter()
        .map(|i| indices_from_coordinates(hexgrid::geo_coordinates(*i), width, total))
        .collect()
}

use std::collections::hash_map::Entry::{Occupied, Vacant};
use std::collections::{BinaryHeap, HashMap};

use std::hash::Hash;

use petgraph::visit::{EdgeRef, IntoEdges, VisitMap, Visitable};
use petgraph::algo::Measure;

// The following code is adapted from the petgraph crate, which is
// Copyright (c) 2015

// Permission is hereby granted, free of charge, to any
// person obtaining a copy of this software and associated
// documentation files (the "Software"), to deal in the
// Software without restriction, including without
// limitation the rights to use, copy, modify, merge,
// publish, distribute, sublicense, and/or sell copies of
// the Software, and to permit persons to whom the Software
// is furnished to do so, subject to the following
// conditions:

// The above copyright notice and this permission notice
// shall be included in all copies or substantial portions
// of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF
// ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED
// TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
// PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT
// SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
// CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
// OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR
// IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//  DEALINGS IN THE SOFTWARE.

use std::cmp::Ordering;

/// `MinScored<K, T>` holds a score `K` and a scored object `T` in
/// a pair for use with a `BinaryHeap`.
///
/// `MinScored` compares in reverse order by the score, so that we can
/// use `BinaryHeap` as a min-heap to extract the score-value pair with the
/// least score.
///
/// **Note:** `MinScored` implements a total order (`Ord`), so that it is
/// possible to use float types as scores.
#[derive(Copy, Clone, Debug)]
pub struct MinScored<K, T>(pub K, pub T);

impl<K: PartialOrd, T> PartialEq for MinScored<K, T> {
    #[inline]
    fn eq(&self, other: &MinScored<K, T>) -> bool {
        self.cmp(other) == Ordering::Equal
    }
}

impl<K: PartialOrd, T> Eq for MinScored<K, T> {}

impl<K: PartialOrd, T> PartialOrd for MinScored<K, T> {
    #[inline]
    fn partial_cmp(&self, other: &MinScored<K, T>) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl<K: PartialOrd, T> Ord for MinScored<K, T> {
    #[inline]
    fn cmp(&self, other: &MinScored<K, T>) -> Ordering {
        let a = &self.0;
        let b = &other.0;
        if a == b {
            Ordering::Equal
        } else if a < b {
            Ordering::Greater
        } else if a > b {
            Ordering::Less
        } else if a.ne(a) && b.ne(b) {
            // these are the NaN cases
            Ordering::Equal
        } else if a.ne(a) {
            // Order NaN less, so that it is last in the MinScore order
            Ordering::Less
        } else {
            Ordering::Greater
        }
    }
}

/// \[Generic\] Dijkstra's shortest path algorithm.
///
/// Compute the length of the shortest path from `start` to every reachable
/// node.
///
/// The graph should be `Visitable` and implement `IntoEdges`. The function
/// `edge_cost` should return the cost for a particular edge, which is used
/// to compute path costs. Edge costs must be non-negative.
///
/// If `goal` is not `None`, then the algorithm terminates once the `goal` node's
/// cost is calculated.
///
/// Returns a `HashMap` that maps `NodeId` to path cost.
/// # Example
/// ```rust
/// use petgraph::Graph;
/// use petgraph::algo::dijkstra;
/// use petgraph::prelude::*;
/// use std::collections::HashMap;
///
/// let mut graph : Graph<(),(),Directed>= Graph::new();
/// let a = graph.add_node(()); // node with no weight
/// let b = graph.add_node(());
/// let c = graph.add_node(());
/// let d = graph.add_node(());
/// let e = graph.add_node(());
/// let f = graph.add_node(());
/// let g = graph.add_node(());
/// let h = graph.add_node(());
/// // z will be in another connected component
/// let z = graph.add_node(());
///
/// graph.extend_with_edges(&[
///     (a, b),
///     (b, c),
///     (c, d),
///     (d, a),
///     (e, f),
///     (b, e),
///     (f, g),
///     (g, h),
///     (h, e)
/// ]);
/// // a ----> b ----> e ----> f
/// // ^       |       ^       |
/// // |       v       |       v
/// // d <---- c       h <---- g
///
/// let expected_res: HashMap<NodeIndex, usize> = [
///      (a, 3),
///      (b, 0),
///      (c, 1),
///      (d, 2),
///      (e, 1),
///      (f, 2),
///      (g, 3),
///      (h, 4)
///     ].iter().cloned().collect();
/// let res = dijkstra(&graph,b,None, |_| 1);
/// assert_eq!(res, expected_res);
/// // z is not inside res because there is not path from b to z.
/// ```
pub fn dijkstra<G, F, K>(
    graph: G,
    start: G::NodeId,
    mut goal: Vec<G::NodeId>,
    mut edge_cost: F,
) -> HashMap<G::NodeId, K>
where
    G: IntoEdges + Visitable,
    G::NodeId: Eq + Hash,
    F: FnMut(G::EdgeRef) -> K,
    K: Measure + Copy,
{
    let mut visited = graph.visit_map();
    let mut scores = HashMap::new();
    //let mut predecessor = HashMap::new();
    let mut visit_next = BinaryHeap::new();
    let zero_score = K::default();
    scores.insert(start, zero_score);
    visit_next.push(MinScored(zero_score, start));
    while let Some(MinScored(node_score, node)) = visit_next.pop() {
        if visited.is_visited(&node) {
            continue;
        }
        if goal.contains(&node) {
            if goal.len() == 1 {
                break;
            }
            let index = goal.iter().position(|&g| g == node).unwrap();
            goal.remove(index);
        }
        for edge in graph.edges(node) {
            let next = edge.target();
            if visited.is_visited(&next) {
                continue;
            }
            let next_score = node_score + edge_cost(edge);
            match scores.entry(next) {
                Occupied(ent) => {
                    if next_score < *ent.get() {
                        *ent.into_mut() = next_score;
                        visit_next.push(MinScored(next_score, next));
                        //predecessor.insert(next.clone(), node.clone());
                    }
                }
                Vacant(ent) => {
                    ent.insert(next_score);
                    visit_next.push(MinScored(next_score, next));
                    //predecessor.insert(next.clone(), node.clone());
                }
            }
        }
        visited.visit(node);
    }
    scores
}

