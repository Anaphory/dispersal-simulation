use std::f64::consts::PI;

pub use libh3::{geo_to_h3, GeoCoord, H3Index as Index};

pub fn closest_grid_point(longitude: f64, latitude: f64) -> Option<Index> {
    let point = GeoCoord {
        lat: latitude * PI / 180.,
        lon: longitude * PI / 180.,
    };
    geo_to_h3(&point, 5).ok()
}


pub fn nearby_locations(index: Index) -> Vec<Index> {
    libh3::k_ring_distances(index, 5)
        .iter()
        .map(|(i, _d)| *i)
        .collect()
}

pub fn geo_coordinates(index: Index) -> GeoCoord {
    libh3::h3_to_geo(index)
}

use lazy_static::lazy_static;
use std::collections::HashMap;
use std::sync::RwLock;
lazy_static! {
    static ref GD_CACHE: RwLock<HashMap<(Index, Index), i32>> = RwLock::new(HashMap::new());
}

pub fn hex_distance(i1: Index, i2: Index) -> i32 {
    match GD_CACHE.read().unwrap().get(&(i1, i2)) {
        None => {
            let d = libh3::h3_distance(i1, i2).unwrap_or(-1);
            match GD_CACHE.try_write() {
                Ok(mut l) => {
                    l.entry((i1, i2)).or_insert(d);
                }
                Err(_) => {}
            }
            d
        }
        Some(d) => *d,
    }
}

// NOT DIRECTLY HEXGRID-RELATED.
// TODO: Move this to a separate module.
use petgraph::visit::{EdgeRef, IntoEdges, VisitMap, Visitable};
use std::hash::Hash;
use petgraph::algo::Measure;
use std::collections::{BinaryHeap};
use std::collections::hash_map::Entry::{Occupied, Vacant};
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

pub fn bounded_dijkstra<G, F, K>(
    graph: G,
    start: G::NodeId,
    max_distance: K,
    mut edge_cost: F
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
        if max_distance < node_score {
            break;
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
