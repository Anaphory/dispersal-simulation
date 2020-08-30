pub use libh3::{geo_to_h3, GeoCoord, H3Index as Index};

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
            if let Ok(mut l) = GD_CACHE.try_write() {
                l.entry((i1, i2)).or_insert(d);
            }
            d
        }
        Some(d) => *d,
    }
}
