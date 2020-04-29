use std::f64::consts::PI;
use rand::prelude::*;
use rand_distr::Normal;
use rayon::prelude::*;
use std::collections::HashMap;

// Hexgrid
use libh3::{geo_to_h3, GeoCoord, H3Index};

fn closest_grid_point(longitude: f64, latitude: f64) -> Option<H3Index> {
    let point = GeoCoord {lat: latitude * PI / 180.,
                          lon: longitude * PI / 180.};
    geo_to_h3(&point, 5).ok()
}

fn nearby_locations(index: H3Index) -> Vec<H3Index> {
    libh3::k_ring_distances(index, 5).iter().map(|(i, _d)| {*i}).collect()
}

fn geo_coordinates(index: H3Index) -> GeoCoord {
    libh3::h3_to_geo(index)
}

// Util

fn patch_from_coordinates(coordinates: GeoCoord, image_pixels: &Vec<u16>, pixels_width: usize) -> Option<KCal> {
    let column = ((coordinates.lon + PI) / PI / 2. * pixels_width as f64).round() as usize;
    let row = ((-coordinates.lat + PI / 2.) / PI * (image_pixels.len() / pixels_width) as f64).round() as usize;
    let index = row * pixels_width + column;
    let precipitation = image_pixels.get(index)?;

    if *precipitation == 0 {
        None
    } else {
        let alpha = (10.0_f32).powf(-8.07);
        // FIXME: 4 is an arbitrary factor
        Some(4. * alpha * (*precipitation as f32).powf(2.64))
    }
}

// Units used
type KCal = f32;
type HalfYears = u32;

struct Family {
    descendence: String,
    culture: Culture,
    location: H3Index,
    location_history: Vec<H3Index>,
    number_offspring: u16,
    effective_size: usize,
    stored_resources: KCal,
    seasons_till_next_child: HalfYears,
    seasons_till_next_mutation: HalfYears
}

impl std::fmt::Debug for Patch {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Patch")
            .field("r", &self.resources)
            .field("/", &self.max_resources)
            .finish()
    }
}

impl std::fmt::Debug for Family {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let geo = geo_coordinates(self.location);

        f.debug_struct("Family")
            .field("descendence", &self.descendence)
            .field("size", &self.effective_size)
            .field("seasons_till_next_child", &self.seasons_till_next_child)
            .field("lat", &(geo.lat * 180. / PI))
            .field("lon", &(geo.lon * 180. / PI))
            .field("culture", &format_args!("{:020b}", &self.culture))
            .field("stored_resources", &self.stored_resources)
            .finish()
    }
}

type Culture = u32;

struct Patch {
    resources: f32,
    max_resources: f32
}

struct State {
    patches: HashMap<H3Index, Patch>,
    //families: HashMap<H3Index, Vec<Family>>,
    families: Vec<Family>,
    t: HalfYears
}

impl PartialEq for Family {fn eq(&self, other: &Self) -> bool {self.descendence == other.descendence}}

fn by_location(families: &Vec<Family>) -> HashMap<H3Index, HashMap<Culture, usize>> {
    let mut cultures_by_location: HashMap<H3Index, HashMap<Culture, usize>> = HashMap::new();
    for family in families.iter() {
        let cultures = cultures_by_location.entry(family.location).or_insert(HashMap::new());
        let counter = cultures.entry(family.culture).or_insert(0);
        *counter += family.effective_size;
    }
    cultures_by_location
}

fn step_part_1(families: &mut Vec<Family>, patches: &HashMap<H3Index, Patch>, p: &Parameters) -> HashMap<H3Index, Vec<Family>> {
    // First, families adjust. These things can happen in parallel for families
    // outside interaction range.
    let mut families_by_location: HashMap<H3Index, Vec<Family>> = HashMap::new();
    let cultures_by_location = by_location(families);

    let mut rng = rand::thread_rng();
    families.shuffle(&mut rng);

    for mut family in families.drain(..) {
        // Submodule 7.1

        use_resources_and_maybe_shrink(&mut family, p);
        if is_moribund(&family) {
            continue;
        }

        // Submodule 7.2
        maybe_grow(&mut family);

        let nearby = nearby_locations(&family.location);
        // Submodule 7.3
        let descendant = maybe_procreate(&mut family);
        match descendant {
            None => {},
            Some(descendant) => {
            // Migration is handled in Submodule 7.4. In terms of scheduling a
            // new family can (and if possible should) move immediately,
            // immediately before their parent. This behaviour is taken from del
            // Castillo (2013). Also, a descendant family will move *before*
            // their progenitor.
                let observed = observe_neighbors(
                    &descendant, patches, &cultures_by_location,
                    &nearby,
                    &p.attention_probability);
                let destination = decide_on_moving(
                    &descendant,
                    observed,
                    true,
                    p);
                match destination {
                    None => {
                        families_by_location.entry(family.location).or_insert(vec![]).push(descendant);
                    }
                    Some(destination) => {
                        families_by_location.entry(destination).or_insert(vec![]).push(descendant);
                    }
                }
            }
        }
        // Subsequent earlier moves affect the possible targets of later moves,
        // so scheduling matters and shuffling is important to remove
        // first-mover effects.

        let observed = observe_neighbors(
            &family, patches, &cultures_by_location,
            &nearby,
            &p.attention_probability);
        let destination = decide_on_moving(
            &family,
            observed,
            true,
            p);
        match destination {
            None => {
                families_by_location.entry(family.location).or_insert(vec![]).push(family);
            },
            Some(destination) => {
                families_by_location.entry(destination).or_insert(vec![]).push(family);
            }
        }
    }
    families_by_location
}

fn step_part_2(families_by_location: &mut HashMap<H3Index, Vec<Family>>, patches: &mut HashMap<H3Index, Patch>, p: &Parameters) {
    // Then, the resources of a patch are updated according to the families
    // exploiting them over the season. This is described in Submodule 7.5.
    // Everything here happens locally to a patch, with no external interaction,
    // so this can be done in parallel.
    for (patch_id, families) in families_by_location {
        let mut patch = patches.entry(*patch_id).or_insert(Patch {resources: 0.0, max_resources: 0.0});
        let mut resource_reduction: KCal = 0.0;
        let (cooperatives, sum_labor) = cooperate(families.iter_mut().collect());

        for cooperating_families in cooperatives {
            resource_reduction += extract_resources(
                &mut patch, cooperating_families, sum_labor,
                p);
        }

        exploit(&mut patch, resource_reduction);
    }

    // Then, patches advance to the next season according to Submodule 7.6.
    for mut patch in patches.values_mut() {
        recover(&mut patch, &p.resource_recovery);
    }

    // The old state has thus been transformed into the new state.
}

fn similar_culture(c1: &Culture, c2: &Culture) -> bool {
    // FIXME wasn't there a parameter for this?
    (c1 ^ c2).count_ones() < 6
}

fn cooperate<'a>(families_in_this_location: Vec<&'a mut Family>) -> (Vec<Vec<&'a mut Family>>, usize){
    let mut sum_labor: usize = 0;
    let mut groups: Vec<Vec<&mut Family>> = vec![];
    for family in families_in_this_location {
        sum_labor += family.effective_size;
        let mut joined_group: Option<&mut Vec<&mut Family>> = None;
        for group in groups.iter_mut() {
            let mut join = true;
            for other_family in group.iter() {
                if !similar_culture(&family.culture, &other_family.culture) {
                    join = false;
                    break;
                }
            }
            if join {
                joined_group = Some(group);
                break;
            }
        }
        match joined_group {
            None => {
                groups.push(vec![family]);
            },
            Some(group) => {
                group.push(family);
            }
        }
    }
    (groups, sum_labor)
}

struct Parameters {
    attention_probability: f32,
    time_step_energy_use: KCal,
    storage_loss: f32,
    resource_recovery: f32,
    culture_mutation_rate: f64,
    culture_dimensionality: u8,
    cooperation_gain: f32,
    accessible_resources: f32,
    evidence_needed: f32,

    boundary_west: f64,
    boundary_east: f64,
    boundary_south: f64,
    boundary_north: f64
}

fn known_location<'a>(history: &Vec<H3Index>, nearby: &Vec<H3Index>, attention_probability: &f32) -> Vec<H3Index> {
    let mut result = vec![];
    for location in nearby {
        if history[0..min(8, history.len())].contains(&location) {
            result.push(*location);
            continue;
        }
        if attention(attention_probability) {
            result.push(*location);
            continue;
        }
    }
    result
}

fn attention(attention_probability: &f32) -> bool {
    random::<f32>() < *attention_probability
}

fn scout<'a>(location: H3Index,
             reference_culture: &Culture,
             patches: &'a HashMap<H3Index, Patch>,
             cultures_by_location: &HashMap<H3Index, HashMap<Culture, usize>>)
             -> Option<(H3Index, &'a Patch, usize, usize)> {
    let cultures = match cultures_by_location.get(&location) {
        None => return Some((location, patches.get(&location)?, 0, 0)),
        Some(c) => c
    };
    let mut cooper: usize = 0;
    let mut compet: usize = 0;
    for (culture, count) in cultures {
        if similar_culture(reference_culture, culture) {
            cooper += count;
        } else {
            compet += count;
        }
    }
    Some((location, patches.get(&location)?, cooper, compet))
}

fn observe_neighbors<'a>(
    family: &'a Family,
    patches: &'a HashMap<H3Index, Patch>,
    cultures_by_location: &HashMap<H3Index, HashMap<Culture, usize>>,
    neighbors: &Vec<H3Index>,
    attention_probability: &f32) -> Vec<(H3Index, &'a Patch, usize, usize)> {
    let mut result: Vec<(H3Index, &'a Patch, usize, usize)> = vec![];
    match scout(family.location, &family.culture, patches, cultures_by_location) {
        None => {},
        Some((location, patch, cooperators, competitors)) =>
            result.push((location, patch, cooperators, competitors))
    }
    for location in known_location(&family.location_history, neighbors, attention_probability) {
        if location == family.location {
            continue
        }
        match scout(location, &family.culture, patches, cultures_by_location) {
            None => {},
            Some((location, patch, cooperators, competitors)) => {
                result.push((location, patch, cooperators, competitors));
            }
        }
    }
    result
}

fn initialization(precipitation: &Vec<u16>, width: usize, p: &Parameters) -> Option<State> {
    let start1: H3Index = closest_grid_point(-159.873 , 65.613)?;
    let start2: H3Index = closest_grid_point(-158.2718, 60.8071)?;

    let area: f32 = libh3::hex_area_km_2(5) as f32;

    let mut patches: HashMap<H3Index, Option<Patch>> = HashMap::new();
    let mut new_patches: std::collections::BinaryHeap<H3Index> = std::collections::BinaryHeap::new();
    new_patches.push(start1);
    new_patches.push(start2);
    loop {
        let next = match new_patches.pop() {
            None => break,
            Some(n) => n
        };
        if patches.contains_key(&next) {
            continue;
        }
        let geo = geo_coordinates(next);
        let latitude = geo.lat * 180. / PI;
        let longitude = geo.lon * 180. / PI;
        if p.boundary_west < longitude && longitude < p.boundary_east &&
            p.boundary_south < latitude && latitude < p.boundary_north {
                let patch = patch_from_coordinates(geo, precipitation, width);
                match patch {
                    None => {
                        patches.insert(next, None);
                        continue;
                    },
                    Some(resources) => {
                        patches.insert(next, Some(Patch{
                            resources: resources * p.time_step_energy_use * area,
                            max_resources: resources * p.time_step_energy_use * area}));
                    }
                };
                for q in nearby_locations(&next) {
                    if patches.contains_key(&q) {
                        continue;
                    }
                    new_patches.push(q);
                }
            } else {
                patches.insert(next, None);
            }
    }
    let f1 = Family {
        descendence: String::from("A"),
        location: start1,
        location_history: vec![],
        seasons_till_next_child: 4,
        culture: 0b000_000_000_000_000,

        effective_size: 2,
        number_offspring: 0,
        seasons_till_next_mutation: 0, // FIXME: Don't mutate immediately
        stored_resources: 1_600_000.
    };
    let f2 = Family {
        descendence: String::from("F"),
        location: start2,
        location_history: vec![],
        seasons_till_next_child: 4,
        culture: 0b111_111_111_111_111,

        effective_size: 2,
        number_offspring: 0,
        seasons_till_next_mutation: 0, // FIXME: Don't mutate immediately
        stored_resources: 1_600_000.
    };
    // let families: HashMap<H3Index, Vec<Family>> = HashMap::new();
    // families.insert(start1, vec![f1]);
    // families.insert(start2, vec![f2]);
    Some(State {
        patches: patches.drain().filter_map(|(i, p)| match p {None => None, Some(q) => Some((i, q))}).collect(),
        families: vec![f1, f2],
        t: 0
    })
}

fn use_resources_and_maybe_shrink(family: &mut Family, p: &Parameters) {
    let resources = family.stored_resources;
    let mut size = family.effective_size;
    while resources_at_season_end(resources, size, p) < 0. && size > 0 {
        size -= 1;
        family.seasons_till_next_child = max(family.seasons_till_next_child, 2)
    }
    family.effective_size = size;
    family.stored_resources = resources_at_season_end(family.stored_resources, size, p);
}

fn resources_at_season_end(resources: KCal, size: usize, p: &Parameters) -> KCal {
    let mut resources_after: KCal = resources -
        (size as f32) * p.time_step_energy_use;
    if resources_after > 0. {
        resources_after *= 1. - p.storage_loss;
    }
    resources_after
}

fn is_moribund(family: &Family) -> bool {
    family.effective_size < 2
}

fn maybe_grow(family: &mut Family) {
    if family.seasons_till_next_child == 0 {
        family.effective_size += 1;
        family.seasons_till_next_child = 2;
    } else {
        family.seasons_till_next_child -= 1;
    }
}

fn maybe_procreate(family: &mut Family) -> Option<Family> {
    if family.effective_size < 10 {
        None
    } else {
        family.number_offspring += 1;
        family.effective_size -= 2;
        Some(Family {
            descendence: format!("{}:{:}", family.descendence, family.number_offspring),
            location: family.location,
            location_history: vec![],
            seasons_till_next_child: 12*2,
            culture: family.culture,

            effective_size: 2,
            number_offspring: 0,
            seasons_till_next_mutation: 0, // FIXME: Don't mutate immediately
            stored_resources: 0.
        })
    }
}

#[test]
fn test_decide_on_moving() {
    let patch1 = Patch{resources: 10., max_resources: 10.};
    let patch2 = Patch{resources: 100., max_resources: 100.};
    let _patch3 = Patch{resources: 100000., max_resources: 100000.};
    let mini_param = Parameters {
            attention_probability: 0.1,
            time_step_energy_use: 10.,
            storage_loss: 0.25,
            resource_recovery: 0.20,
            culture_mutation_rate: 6e-3,
            culture_dimensionality: 20,
            cooperation_gain: 0.0,
            accessible_resources: 1.0,
            evidence_needed: 0.3,

        // Some part of Western Alaska
            boundary_west: -168.571541,
            boundary_east: -148.571541,
            boundary_south: 56.028198,
            boundary_north: 74.52671,
    };

    let mini_family = Family {
        culture: 0,
        descendence: String::from(""),
        effective_size: 2,
        location: 0,
        location_history: vec![1],
        number_offspring: 0,
        seasons_till_next_child: 0,
        seasons_till_next_mutation: 0,
        stored_resources: 1.
    };
    assert_eq!(decide_on_moving(
        &mini_family,
        vec![(1, &patch1, mini_family.effective_size, 0)],
        false,
        &mini_param
    ), Some(1));

    assert_eq!(decide_on_moving(
        &mini_family,
        vec![(1, &patch1, mini_family.effective_size, 0),
             (2, &patch2, mini_family.effective_size, 0)],
        true,
        &mini_param
    ), Some(2));

    assert_eq!(decide_on_moving(
        &mini_family,
        vec![(1, &patch1, mini_family.effective_size, 0),
             (2, &patch2, mini_family.effective_size, 0)],
        false,
        &mini_param
    ), Some(2));
}

#[test]
fn test_decide_on_moving_is_uniform() {
    let patch1 = Patch{resources: 10., max_resources: 10.};
    let patch2 = Patch{resources: 100., max_resources: 100.};
    let patch3 = Patch{resources: 100000., max_resources: 100000.};
    let mini_param = Parameters {
            attention_probability: 0.1,
            time_step_energy_use: 10.,
            storage_loss: 0.25,
            resource_recovery: 0.20,
            culture_mutation_rate: 6e-3,
            culture_dimensionality: 20,
            cooperation_gain: 0.0,
            accessible_resources: 1.0,
            evidence_needed: 0.3,

        // Some part of Western Alaska
            boundary_west: -168.571541,
            boundary_east: -148.571541,
            boundary_south: 56.028198,
            boundary_north: 74.52671,
    };

    let mini_family = Family {
        culture: 0,
        descendence: String::from(""),
        effective_size: 2,
        location: 0,
        location_history: vec![1],
        number_offspring: 0,
        seasons_till_next_child: 0,
        seasons_till_next_mutation: 0,
        stored_resources: 1.
    };

    let mut c = [0, 0, 0, 0, 0];
    for _ in 1..300 {
        let k = decide_on_moving(
            &mini_family,
            vec![(1, &patch1, mini_family.effective_size, 0),
                 (2, &patch3, mini_family.effective_size, 0),
                 (3, &patch2, mini_family.effective_size, 0),
                 (4, &patch2, mini_family.effective_size, 0)],
            false,
            &mini_param
        );
        c[match k {None => 0, Some(k) => k as usize}] += 1;
    }
    assert_eq!(c[0], 0);
    assert_eq!(c[1], 0);
    assert!((c[2] - 100i8).abs() < 15);
    assert!((c[3] - 100i8).abs() < 15);
    assert!((c[4] - 100i8).abs() < 15);
}

fn decide_on_moving<'a>(
    family: &'a Family,
    known_destinations: Vec<(H3Index, &Patch, usize, usize)>,
    avoid_stay: bool,
    p: &Parameters) -> Option<H3Index> {
    let mut kd = known_destinations.iter();
    let (mut target, mut patch, mut cooperators, mut competitors): (H3Index, &Patch, usize, usize);
    match kd.next() {
        None => {
            return None
        }
        Some((t, p, cooper, compet)) => {
            target = *t;
            patch = p;
            cooperators = *cooper;
            competitors = *compet;
        }
    };
    if avoid_stay {
        match kd.next() {
            None => {
                return None
            }
            Some((t, p, cooper, compet)) => {
                target = *t;
                patch = p;
                cooperators = *cooper;
                competitors = *compet;
            }
        }
    }
    let mut max_gain: KCal = resources_from_patch(patch, cooperators, competitors, true, p) *
        family.effective_size as f32 / cooperators as f32;
    let current_gain = if avoid_stay {0.} else {max_gain};
    let mut c = 0.;
    for (coords, patch, cooperators, competitors) in kd {
        let expected_gain = resources_from_patch(
            patch,
            family.effective_size + *cooperators,
            *competitors,
            true,
            p
        ) * (family.effective_size as f32/ (family.effective_size + *cooperators) as f32);
        if expected_gain >= max_gain {
            if expected_gain <
                current_gain +
                p.time_step_energy_use * p.evidence_needed {
                    continue;
                }
            if (expected_gain - max_gain).abs() < std::f32::EPSILON {
                c += 1.;
                if random::<f32>() < c / (1. + c) {
                    continue
                }
            }
            target = *coords;
            max_gain = expected_gain;
        }
    }
    Some(target)
}

fn extract_resources(patch: &mut Patch, group: Vec<&mut Family>, total_labor_here: usize, p: &Parameters) -> KCal {
    let labor: usize = group.iter().map(|f| {f.effective_size as usize}).sum();
    let resources_extracted = resources_from_patch(
        &patch, labor, total_labor_here - labor,
        false,
        &p);
    let mut group_copy: Vec<&mut Family> = vec![];
    for family in group {
        family.stored_resources +=
            resources_extracted * family.effective_size as f32 / labor as f32;
        group_copy.push(family);
    }
    adjust_culture(group_copy, p);
    resources_extracted
}

#[test]
fn test_resources_from_patch() {
    let patch1 = Patch{resources: 10., max_resources: 10.};
    let patch2 = Patch{resources: 100., max_resources: 100.};
    let patch3 = Patch{resources: 100000., max_resources: 100000.};
    let mut mini_param = Parameters {
            attention_probability: 0.1,
            time_step_energy_use: 10.,
            storage_loss: 0.25,
            resource_recovery: 0.20,
            culture_mutation_rate: 6e-3,
            culture_dimensionality: 20,
            cooperation_gain: 0.0,
            accessible_resources: 1.0,
            evidence_needed: 0.3,

        // Some part of Western Alaska
            boundary_west: -168.571541,
            boundary_east: -148.571541,
            boundary_south: 56.028198,
            boundary_north: 74.52671,
    };
    assert_eq!(resources_from_patch(&patch1, 2, 0, true, &mini_param), 10.);
    assert_eq!(resources_from_patch(&patch2, 2, 0, true, &mini_param), 22.);
    assert_eq!(resources_from_patch(&patch3, 2, 0, true, &mini_param), 22.);
    assert!((resources_from_patch(&patch2, 3, 7, true, &mini_param) - 30.).abs() < 0.001);
    mini_param.cooperation_gain = 0.5;
    assert!((resources_from_patch(&patch3, 25, 0, true, &mini_param) - 1.5*25.*10.).abs() < 5.);
}

fn resources_from_patch(patch: &Patch, labor: usize, others_labor: usize, estimate: bool, p: &Parameters) -> KCal {
    let mut my_relative_returns: f32 =
        p.time_step_energy_use * labor as f32 *
        effective_labor_through_cooperation(labor, &p.cooperation_gain);
    let mut rng = rand::thread_rng();
    if !estimate {
        let dist = Normal::new(my_relative_returns, p.time_step_energy_use / (labor as f32).powf(0.5));
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
            p.time_step_energy_use * others_labor as f32 *
            effective_labor_through_cooperation(labor, &p.cooperation_gain);
        if !estimate {
            let dist = Normal::new(my_relative_returns, p.time_step_energy_use / (labor as f32).powf(0.5));
            match dist {
                Err(_) => {}
                Ok(d) => {
                    others_relative_returns = d.sample(&mut rng);
                    others_relative_returns = max(0., others_relative_returns)
                }
            }
        }
    }
    my_relative_returns /
        (my_relative_returns + others_relative_returns) * min(
            my_relative_returns + others_relative_returns,
            patch.resources * p.accessible_resources
        )
}

fn effective_labor_through_cooperation(n_cooperators: usize, cooperation_gain: &f32) -> f32 {
    1. + (n_cooperators as f32 - 1.).powf(*cooperation_gain) / 10.
}

fn adjust_culture(mut cooperating_families: Vec<&mut Family>, p: &Parameters) {
    for f in cooperating_families.iter_mut() {
        mutate_culture(f, p);
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

fn mutate_culture(family: &mut Family, p: &Parameters) {
    if family.seasons_till_next_mutation > 0 {
        family.seasons_till_next_mutation -= 1;
    }
    if family.seasons_till_next_mutation == 0 {
        let i: u8 = rand::thread_rng().gen_range(0, p.culture_dimensionality);
        family.culture ^= 1 << i;
        family.seasons_till_next_mutation = random::<f64>().log(1. - p.culture_mutation_rate) as u32;
    }
}

fn exploit(patch: &mut Patch, resource_reduction: KCal) {
    patch.resources -= resource_reduction;
    //FIXME why does this keep happening?
    // assert!(patch.resources > 0.);
}

fn recover(patch: &mut Patch, resource_recovery: &f32) {
    if patch.resources < patch.max_resources - 1. {
        patch.resources += patch.resources *
            resource_recovery *
            (1. - patch.resources / patch.max_resources);
    }
}

fn min<T: PartialOrd>(a: T, b: T) -> T {if a < b {a} else {b}}
fn max<T: PartialOrd>(a: T, b: T) -> T {if b < a {a} else {b}}

fn main() {
    run();
}

fn run() -> Option<()> {
    let mut parser = argparse::ArgumentParser::new();
    parser.set_description("Run a dispersal simulation");
    let p = Parameters {
        attention_probability: 0.1,
        time_step_energy_use: 2263. as KCal * 365.242_2 / 2.,
        storage_loss: 0.33,
        resource_recovery: 0.20,
        culture_mutation_rate: 6e-3,
        culture_dimensionality: 20,
        cooperation_gain: 0.5,
        accessible_resources: 0.2,
        evidence_needed: 0.3,

        boundary_west: -168.571541,
        boundary_east: -34.535395,
        boundary_south: -56.028198,
        boundary_north: 74.52671,
    };
    let f = std::fs::File::open("wc2.1_5m_bio_12-16bit.tif").ok()?;
    let maybe_image = tiff::decoder::Decoder::new(f);
    let mut image = maybe_image.ok()?;
    let (width, height) = image.dimensions().ok()?;
    println!("# {}x{}", width, height);
    let outcome = image.read_image().ok()?;
    println!("# Image read");
    let vec = match outcome {
        tiff::decoder::DecodingResult::U8(v) => v.iter().map(|g| {*g as u16}).collect(),
        tiff::decoder::DecodingResult::U16(w) => w
    };
    println!("# Initialization…");
    // FIXME: There is something messed up in lat/long -> image pixels. It works
    // currently, but the names point towards a deep issue with the layout of
    // the tiff in memory or similar.
    let mut s: State = initialization(&vec, width as usize, &p)?;

    loop {
        let mut families_by_location = step_part_1(&mut s.families, &s.patches, &p);
        println!("{:?}", families_by_location.iter().map(|(k, v)| {
            let g = geo_coordinates(*k);
            (g.lon, g.lat, v.len())
        }).collect::<Vec<(f64, f64, usize)>>());
        step_part_2(&mut families_by_location, &mut s.patches, &p);
        for (location, families) in families_by_location{
            for mut family in families {
                family.location_history.push(family.location);
                family.location = location;
                s.families.push(family);
            }
        };
        println!("{:}", s.families.len());
        if s.families.is_empty() {
            break;
        }
        s.t += 1;
        println!("t={:}", s.t);
    }
    Some(())
}
