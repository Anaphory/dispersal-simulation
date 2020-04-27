use std::f64::consts::PI;
use rand::prelude::*;
use rand_distr::Normal;
use rayon::prelude::*;

// Hexgrid
use libh3::{geo_to_h3, GeoCoord, H3Index};

fn closest_grid_point(longitude: f64, latitude: f64) -> Option<H3Index> {
    let point = GeoCoord {lat: latitude* PI / 180.,
                          lon: longitude * PI / 180.};
    return geo_to_h3(&point, 5).ok();
}

fn geo_coordinates(index: H3Index) -> GeoCoord {
    return libh3::h3_to_geo(index);
}

// Util

fn patch_from_coordinates(coordinates: GeoCoord, data: &Vec<u16>, pixels_width: usize) -> Option<KCal> {
    let column = ((coordinates.lon + PI) / PI / 2. * pixels_width as f64).round() as usize;
    let row = ((coordinates.lat + PI / 2.) / PI * (data.len() / pixels_width) as f64).round() as usize;
    let index = row * data.len() / pixels_width + column;
    let precipitation = data.get(index)?;
    if *precipitation == 0 {
        return None;
    } else {
        let alpha = (10.0_f32).powf(-8.07);
        return Some(alpha * (*precipitation as f32).powf(2.64));
    }
}

// Units used
type KCal = f32;
type HalfYears = u32;

static attention_probability: f32 = 0.1;
static time_step_energy_use: KCal = 10.;
static storage_loss: f32 = 0.33;
static resource_recovery: f32 = 0.20;
static culture_mutation_rate: f32 = 6e-3;
static culture_dimensionality: u8 = 20;
static cooperation_gain: f32 = 0.5;
static accessible_resources: f32 = 0.2;
static evidence_needed: f32 = 0.3;

static boundary_west: f64 = -168.571541;
static boundary_east: f64 = -34.535395;
static boundary_south: f64 = -56.028198;
static boundary_north: f64 = 74.52671;

struct Family {
    descendence: String,
    culture: Culture,
    location_history: Vec<H3Index>,
    number_offspring: u16,
    effective_size: usize,
    stored_resources: KCal,
    seasons_till_next_child: HalfYears,
    seasons_till_next_mutation: HalfYears
}

type Culture = u32;

struct Patch {
    resources: f32,
    max_resources: f32
}

struct State {
    patches: std::collections::HashMap<H3Index, Patch>,
    families: Vec<Family>,
    t: HalfYears
}



fn step(state: &mut State) {
    let mut rng = thread_rng();

    state.families.shuffle(&mut rng);
    for family in &mut state.families {
        use_resources_and_maybe_shrink(family);
    }
}

fn known_location<'a>(history: &'a Vec<H3Index>, nearby: &'a Vec<H3Index>) -> Vec<&'a H3Index> {
    return history[0..4].iter().chain(
        nearby.iter().filter(|n| {
            !history[0..4].contains(n) && (
                history[4..8].contains(n) || attention())}
        )
    ).collect();
}

fn attention() -> bool {
    return random::<f32>() < attention_probability;
}

fn observe_neghbors<'a>(
    family: &'a Family,
    patches: &'a std::collections::HashMap<H3Index, Patch>,
    all_families: std::collections::HashMap<H3Index, Vec<Family>>,
    neighbors: &'a Vec<H3Index>) -> Vec<(&'a H3Index, Option<&'a Patch>, usize, usize)> {
    known_location(&family.location_history, neighbors).iter().filter(
        |k| {patches.contains_key(k)}
    ).map(|l| {
        // FIXME incomplete
        (*l, patches.get(l), 0, 0)
    }).collect()
}

fn initialization(data: &Vec<u16>, width: usize) -> Option<State> {
    let start1: H3Index = closest_grid_point(-159.873 , 65.613)?;
    let start2: H3Index = closest_grid_point(-158.2718, 60.8071)?;

    let mut patches: std::collections::HashMap<H3Index, Patch> = std::collections::HashMap::new();
    let mut new_patches: std::collections::HashSet<H3Index> = std::collections::HashSet::new();
    new_patches.insert(start1);
    new_patches.insert(start2);
    while !new_patches.is_empty() {
        let next = match new_patches.drain().next() {
            None => break,
            Some(n) => n
        };
        if patches.contains_key(&next) {
            continue;
        }
        let geo = geo_coordinates(next);
        let latitude = geo.lat * 180. / PI;
        let longitude = geo.lon * 180. / PI;
        println!("{:}: {:}, {:}", next, longitude, latitude);
        if boundary_west < longitude && longitude < boundary_east &&
            boundary_south < latitude && latitude < boundary_north {
                let p = patch_from_coordinates(geo, data, width);
                match p {
                    None => continue,
                    Some(resources) => patches.insert(next, Patch{
                        resources: resources,
                        max_resources: resources,
                    })
                };
                for (q, r) in libh3::k_ring_distances(next, 5) {
                    if patches.contains_key(&q) {
                        continue;
                    }
                    new_patches.push(q);
                }
            }
    }
    return Some(State {
        patches: patches,
        families: vec![
            Family {
                descendence: String::from("A"),
                location_history: vec![start1],
                seasons_till_next_child: 12*2,
                culture: 0b000_000_000_000_000,

                effective_size: 2,
                number_offspring: 0,
                seasons_till_next_mutation: 0, // FIXME: Don't mutate immediately
                stored_resources: 16_000_000.
            },
            Family {
                descendence: String::from("F"),
                location_history: vec![start2],
                seasons_till_next_child: 12*2,
                culture: 0b111_111_111_111_111,

                effective_size: 2,
                number_offspring: 0,
                seasons_till_next_mutation: 0, // FIXME: Don't mutate immediately
                stored_resources: 16_000_000.
            }
        ],
        t: 0
    })
}

fn use_resources_and_maybe_shrink(family: &mut Family) {
    let resources = family.stored_resources;
    let mut size = family.effective_size;
    while resources_at_season_end(resources, size) < 0. && size > 0 {
        size -= 1;
        family.seasons_till_next_child = std::cmp::max(family.seasons_till_next_child, 2)
    }
    family.effective_size = size;
    family.stored_resources = resources_at_season_end(family.stored_resources, size);
}

fn resources_at_season_end(resources: KCal, size: usize) -> KCal {
    let mut resources_after: KCal = resources -
        (size as f32) * time_step_energy_use;
    if resources_after > 0. {
        resources_after *= 1. - storage_loss;
    }
    return resources_after;
}

fn is_moribund(family: Family) -> bool {
    if family.effective_size < 2 {
        return true;
    } else {
        return false;
    }
}

fn maybe_grow(family: &mut Family) {
    if family.seasons_till_next_child <= 0 {
        family.effective_size += 1;
        family.seasons_till_next_child = 2;
    } else {
        family.seasons_till_next_child -= 1;
    }
}

fn maybe_procreate(family: &mut Family) -> Option<Family> {
    if family.effective_size < 10 {
        return None
    } else {
        family.number_offspring += 1;
        family.effective_size -= 2;
        return Some(Family {
            descendence: String::from("string"),
            location_history: vec![*family.location_history.get(0)?],
            seasons_till_next_child: 12*2,
            culture: family.culture,

            effective_size: 2,
            number_offspring: 0,
            seasons_till_next_mutation: 0, // FIXME: Don't mutate immediately
            stored_resources: 0.
        });
    }
}

fn decide_on_moving<'a>(
    family: &'a Family,
    current_patch: Patch,
    known_destinations: Vec<(H3Index, Patch, usize, usize)>,
    avoid_stay: bool) -> Option<H3Index> {
    let mut kd = known_destinations.iter();
    let (mut target, mut patch, mut cooperators, mut competitors): (H3Index, &Patch, usize, usize);
    match kd.next() {
        None => {
            return None
        }
        Some((t, p, coop, comp)) => {
            target = *t;
            patch = p;
            cooperators = *coop;
            competitors = *comp;
        }
    };
    if avoid_stay {
        match kd.next() {
            None => {
                return None
            }
            Some((t, p, coop, comp)) => {
                target = *t;
                patch = p;
                cooperators = *coop;
                competitors = *comp;
            }
        }
    }
    let mut max_gain: KCal = resources_from_patch(patch, cooperators, competitors, true) *
        family.effective_size as f32 / cooperators as f32;
    let current_gain = if avoid_stay {0.} else {max_gain};
    let mut c = 0.;
    for (coords, patch, cooperators, competitors) in kd {
        let expected_gain = resources_from_patch(
            patch,
            family.effective_size + *cooperators,
            *competitors,
            true
        ) * (family.effective_size / (family.effective_size + *cooperators)) as f32;
        if expected_gain >= max_gain {
            if expected_gain <
                current_gain +
                time_step_energy_use * evidence_needed {
                    continue;
                }
            if expected_gain == max_gain {
                c += 1.;
                if random::<f32>() < c / (1. + c) {
                    continue
                }
            }
            target = *coords;
            max_gain = expected_gain;
        }
    }
    return Some(target);
}


fn movement_bookkeeping(family: &mut Family, destination: H3Index, state: &mut State) -> Option<()> {
    if destination == *family.location_history.get(0)? {
        return Some(());
    }
    family.location_history.push(destination);
    return Some(());
}


fn extract_resources(patch: Patch, group: Vec<&mut Family>, total_labor_here: usize) -> KCal {
    let labor: u8 = group.iter().map(|f| {f.effective_size as u8}).sum();
    let resources_extracted = resources_from_patch(
        &patch, labor as usize, total_labor_here - labor as usize,
        false);
    for family in group {
        family.stored_resources +=
            resources_extracted * family.effective_size as f32 / labor as f32
    }
    return resources_extracted
}

fn resources_from_patch(patch: &Patch, labor: usize, others_labor: usize, estimate: bool) -> KCal {
    let mut my_relative_returns: f32 =
        time_step_energy_use * labor as f32 *
        effective_labor_through_cooperation(labor);
    let mut rng = rand::thread_rng();
    if !estimate {
        let dist = Normal::new(my_relative_returns, time_step_energy_use / (labor as f32).powf(0.5));
        match dist {
            Err(_) => {}
            Ok(d) => {
                my_relative_returns = d.sample(&mut rng);
                my_relative_returns = max(0., my_relative_returns)
            }
        }
    }
    let mut others_relative_returns: f32 = 0.;
    if others_labor > 0 {
        others_relative_returns =
            time_step_energy_use * others_labor as f32 *
            effective_labor_through_cooperation(labor);
        if !estimate {
            let dist = Normal::new(my_relative_returns, time_step_energy_use / (labor as f32).powf(0.5));
            match dist {
                Err(_) => {}
                Ok(d) => {
                    others_relative_returns = d.sample(&mut rng);
                    others_relative_returns = max(0., others_relative_returns)
                }
            }
        }
    }
    return my_relative_returns /
        (my_relative_returns + others_relative_returns) * min(
            my_relative_returns + others_relative_returns,
            patch.resources * accessible_resources
        );
}


fn effective_labor_through_cooperation(n_cooperators: usize) -> f32 {
    return 1. + (n_cooperators as f32 - 1.).powf(cooperation_gain) / 10.
}


fn adjust_culture(mut cooperating_families: Vec<&mut Family>) {
    for f in cooperating_families.iter_mut() {
        mutate_culture(f);
    }
    let family: Option<&&mut Family> = cooperating_families.choose(&mut rand::thread_rng());
    match family {
        None => {},
        Some(f) => {
            let target: Culture = f.culture;
            for f2 in cooperating_families {
                f2.culture = target;
            }
        }
    }
}

fn mutate_culture(family: &mut Family) {
    if family.seasons_till_next_mutation > 0 {
        family.seasons_till_next_mutation -= 1;
    }
    if family.seasons_till_next_mutation == 0 {
        let i: u8 = rand::thread_rng().gen_range(0, culture_dimensionality);
        family.culture ^= 1 << i;
    }
}

fn exploit(patch: &mut Patch, resource_reduction: KCal) {
    patch.resources -= resource_reduction;
    assert!(patch.resources > 0.);
}

fn recover(patch: &mut Patch) {
    if patch.resources < patch.max_resources - 1. {
        patch.resources += patch.resources *
            resource_recovery *
            (1. - patch.resources / patch.max_resources);
    }
}

fn min(a: f32, b: f32) -> f32 {if a < b {a} else {b}}
fn max(a: f32, b: f32) -> f32 {if a > b {a} else {b}}

fn main() {
    run();
}

fn run() -> Option<()> {
    let f = std::fs::File::open("wc2.1_5m_bio_12-16bit.tif").ok()?;
    let maybe_image = tiff::decoder::Decoder::new(f);
    let mut image = maybe_image.ok()?;
    let (width, height) = image.dimensions().ok()?;
    println!("{}x{}", width, height);
    let outcome = image.read_image().ok()?;
    println!("Image read");
    let vec = match outcome {
        tiff::decoder::DecodingResult::U8(v) => v.iter().map(|g| {*g as u16}).collect(),
        tiff::decoder::DecodingResult::U16(w) => w
    };
    println!("Initialization…");
    let s: State = initialization(&vec, width as usize)?;
    println!("Checking patches");

    // println!("x = {:?}", vec);
    for p in s.patches.keys() {
        let g = geo_coordinates(*p);
        println!("[{}, {}],", g.lon, g.lat);
    }

    return Some(());
}


