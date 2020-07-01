use minimal_metabolic_cost_surface::KCal;
use std::path::Path;

pub fn load_precipitation_tif() -> Option<(Vec<u32>, u32)> {
    minimal_metabolic_cost_surface::load_tif(Path::new("wc2.1_5m_bio_12-16bit.tif"))
}


pub fn precipitation_to_resources(
    precipitation: u32
) -> Option<KCal> {
    if precipitation == 0 {
        None
    } else {
        let alpha = (10.0_f32).powf(-8.07);
        let beta: f32 = 2.64;
        // FIXME: 4 is an arbitrary factor
        Some(4. * alpha * (precipitation as f32).powf(beta))
    }
}
