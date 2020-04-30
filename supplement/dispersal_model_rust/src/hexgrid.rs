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
