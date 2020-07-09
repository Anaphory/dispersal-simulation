use crate::hexgrid;
use crate::KCal;
use std::f64::consts::PI;

pub fn load_precipitation_tif() -> Option<(Vec<u16>, u32)> {
    // TODO: How do I specify data paths?
    let f = std::fs::File::open("wc2.1_5m_bio_12-16bit.tif").ok()?;
    let maybe_image = tiff::decoder::Decoder::new(f);
    let mut image = maybe_image.ok()?;
    let (width, height) = image.dimensions().ok()?;
    println!("# {}x{}", width, height);
    let outcome = image.read_image().ok()?;
    println!("# Image read");
    let vec = match outcome {
        tiff::decoder::DecodingResult::U8(v) => v.iter().map(|g| u16::from(*g)).collect(),
        tiff::decoder::DecodingResult::U16(w) => w,
    };
    Some((vec, width))
}

pub fn patch_from_coordinates(
    coordinates: hexgrid::GeoCoord,
    image_pixels: &Vec<u16>,
    pixels_width: usize,
) -> Option<KCal> {
    let column = ((coordinates.lon + PI) / PI / 2. * pixels_width as f64).round() as usize;
    let row = ((-coordinates.lat + PI / 2.) / PI * (image_pixels.len() / pixels_width) as f64)
        .round() as usize;
    let index = row * pixels_width + column;
    let precipitation = image_pixels.get(index)?;

    if *precipitation == 0 {
        None
    } else {
        let alpha = (10.0_f32).powf(-8.07);
        let beta: f32 = 2.64;
        // FIXME: 4 is an arbitrary factor
        Some(4. * alpha * (*precipitation as f32).powf(beta))
    }
}
