use crate::*;

#[test]
pub fn test_decide_on_moving() {
    let patch1 = Patch {
        resources: vec![(0, (10., 10.))].drain(..).collect(),
    };
    let patch2 = Patch {
        resources: vec![(0, (100., 100.))].drain(..).collect(),
    };
    let _patch3 = Patch {
        resources: vec![(0, (10000., 10000.))].drain(..).collect(),
    };
    let mini_param = Parameters {
        attention_probability: 0.1,
        time_step_energy_use: 10.,
        storage_loss: 0.25,
        resource_recovery: 0.20,
        culture_mutation_rate: 6e-3,
        culture_dimensionality: 20,
        cooperation_threshold: 6,
        cooperation_gain: 0.0,
        accessible_resources: 1.0,
        evidence_needed: 0.3,
        payoff_std: 0.1,

        // Some part of Western Alaska
        boundary_west: -168.571541,
        boundary_east: -148.571541,
        boundary_south: 56.028198,
        boundary_north: 74.52671,

        dispersal_graph: petgraph::csr::Csr::new(),
    };

    let mini_family = Family {
        culture: 0,
        descendence: String::from(""),
        effective_size: 2,
        location: 0,
        location_history: vec![1],
        number_offspring: 0,
        seasons_till_next_child: 0,
        seasons_till_next_mutation: None,
        stored_resources: 1.,
        adaptation: ecology::Ecovector::default()
    };
    assert_eq!(
        adaptation::decide_on_moving(
            &mini_family,
            vec![(1, &patch1, mini_family.effective_size, 0)],
            false,
            &mini_param
        ),
        Some(1)
    );

    assert_eq!(
        adaptation::decide_on_moving(
            &mini_family,
            vec![
                (1, &patch1, mini_family.effective_size, 0),
                (2, &patch2, mini_family.effective_size, 0)
            ],
            true,
            &mini_param
        ),
        Some(2)
    );

    assert_eq!(
        adaptation::decide_on_moving(
            &mini_family,
            vec![
                (1, &patch1, mini_family.effective_size, 0),
                (2, &patch2, mini_family.effective_size, 0)
            ],
            false,
            &mini_param
        ),
        Some(2)
    );
}

#[test]
pub fn test_decide_on_moving_is_uniform() {
    let patch1 = Patch {resources: vec![(0, (10., 10.))].drain(..).collect()};
    let patch2 = Patch {resources: vec![(0, (100., 100.))].drain(..).collect()};
    let patch3 = Patch {resources: vec![(0, (10000., 10000.))].drain(..).collect()};
    let mini_param = Parameters {
        attention_probability: 0.1,
        time_step_energy_use: 10.,
        storage_loss: 0.25,
        resource_recovery: 0.20,
        culture_mutation_rate: 6e-3,
        culture_dimensionality: 20,
        cooperation_threshold: 6,
        cooperation_gain: 0.0,
        accessible_resources: 1.0,
        evidence_needed: 0.3,
        payoff_std: 0.1,

        // Some part of Western Alaska
        boundary_west: -168.571541,
        boundary_east: -148.571541,
        boundary_south: 56.028198,
        boundary_north: 74.52671,

        dispersal_graph: petgraph::csr::Csr::new(),
    };

    let mini_family = Family {
        culture: 0,
        descendence: String::from(""),
        effective_size: 2,
        location: 0,
        location_history: vec![1],
        number_offspring: 0,
        seasons_till_next_child: 0,
        seasons_till_next_mutation: None,
        stored_resources: 1.,
        adaptation: ecology::Ecovector::default()
    };

    let mut c = [0, 0, 0, 0, 0];
    for _ in 1..300 {
        let k = adaptation::decide_on_moving(
            &mini_family,
            vec![
                (1, &patch1, mini_family.effective_size, 0),
                (2, &patch3, mini_family.effective_size, 0),
                (3, &patch2, mini_family.effective_size, 0),
                (4, &patch2, mini_family.effective_size, 0),
            ],
            false,
            &mini_param,
        );
        c[match k {
            None => 0,
            Some(k) => k as usize,
        }] += 1;
    }
    assert_eq!(c[0], 0);
    assert_eq!(c[1], 0);
    assert!((c[2] - 100i8).abs() < 15);
    assert!((c[3] - 100i8).abs() < 15);
    assert!((c[4] - 100i8).abs() < 15);
}

