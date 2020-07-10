use petgraph::visit::{EdgeRef, IntoEdges, IntoNeighbors, VisitMap, Visitable, IntoEdgeReferences, Data, GraphRef, GraphBase};
use std::hash::Hash;
use petgraph::algo::Measure;
use std::collections::{BinaryHeap, HashMap, HashSet};
use std::collections::hash_map::Entry::{Occupied, Vacant};
use std::cmp::Ordering;
use rusqlite::{params, Connection, Result, OpenFlags};
use crate::hexgrid;

pub fn edge_costs_from_db(h: &hexgrid::Index) -> Result<Vec<(hexgrid::Index, f64)>> {
    let conn = Connection::open_with_flags(
        "/home/gereon/Public/settlement-of-americas/supplement/distances/plot.sqlite",
        OpenFlags::SQLITE_OPEN_READ_ONLY)?;

    let mut stmt = conn.prepare("SELECT hexbin1, hexbin2, distance FROM dist WHERE hexbin1 = ?")?;
    let result = match stmt.query_map(
        params![*h as i64],
        |row| Ok((row.get::<usize, i64>(1)? as u64, row.get::<usize, f64>(2)?))
    ) {
        Ok(k) => {
            Ok(k.filter_map(Result::ok).collect())
        }
        Err(e) => {
            Err(e)
        }
    };
    result
}

// From petgraph

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

pub fn bounded_dijkstra<N, F, K>(
    start: N,
    max_distance: K,
    edge_cost: F
) -> HashMap<N, K>
where
    N: Eq + Hash + Copy,
    F: Fn(&N) -> Result<Vec<(N, K)>>,
    K: Measure + Copy,
{
    let mut visited: HashSet<N> = HashSet::new();
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
        for (next, cost) in edge_cost(&node).unwrap_or(vec![]) {
            if visited.is_visited(&next) {
                continue;
            }
            let next_score = node_score + cost;
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
