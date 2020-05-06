use crate::*;

#[test]
pub fn test_decide_on_moving() {
    let patch1 = Patch {
        resources: 10.,
        max_resources: 10.,
    };
    let patch2 = Patch {
        resources: 100.,
        max_resources: 100.,
    };
    let _patch3 = Patch {
        resources: 100000.,
        max_resources: 100000.,
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
    let patch1 = Patch {
        resources: 10.,
        max_resources: 10.,
    };
    let patch2 = Patch {
        resources: 100.,
        max_resources: 100.,
    };
    let patch3 = Patch {
        resources: 100000.,
        max_resources: 100000.,
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

#[test]
pub fn test_resources_from_patch() {
    let patch1 = Patch {
        resources: 10.,
        max_resources: 10.,
    };
    let patch2 = Patch {
        resources: 100.,
        max_resources: 100.,
    };
    let patch3 = Patch {
        resources: 100000.,
        max_resources: 100000.,
    };
    let mut mini_param = Parameters {
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
    };
    assert_eq!(
        submodels::ecology::resources_from_patch(&patch1, 2., 0., true, &mini_param),
        10.
    );
    assert_eq!(
        submodels::ecology::resources_from_patch(&patch2, 2., 0., true, &mini_param),
        22.
    );
    assert_eq!(
        submodels::ecology::resources_from_patch(&patch3, 2., 0., true, &mini_param),
        22.
    );
    assert!(
        (submodels::ecology::resources_from_patch(&patch2, 3., 7., true, &mini_param) - 30.).abs()
            < 0.001
    );
    mini_param.cooperation_gain = 0.5;
    assert!(
        (submodels::ecology::resources_from_patch(&patch3, 25., 0., true, &mini_param)
            - 1.5 * 25. * 10.)
            .abs()
            < 5.
    );

    assert!(
        submodels::ecology::resources_from_patch(&patch3, 5., 20., false, &mini_param)
            + submodels::ecology::resources_from_patch(&patch3, 5., 20., false, &mini_param)
            + submodels::ecology::resources_from_patch(&patch3, 5., 20., false, &mini_param)
            + submodels::ecology::resources_from_patch(&patch3, 5., 20., false, &mini_param)
            + submodels::ecology::resources_from_patch(&patch3, 5., 20., false, &mini_param)
            < submodels::ecology::resources_from_patch(&patch3, 25., 0., false, &mini_param)
    )
}
