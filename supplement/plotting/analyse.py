import numpy
import bisect
import pandas
import sklearn
import itertools
import sklearn.ensemble
from matplotlib import pyplot as plt
from osm import areas

from tqdm import tqdm
from interpret_data import data, times, populations

all_predictors = [
    "culture_mutation_rate",
    "culture_dimensionality",
    "cooperation_threshold",
    "minimum_adaptation",
    "fight_deadliness",
    "maximum_resources_one_adult_can_harvest",
    "enemy_discount",
    "season_length_in_years",
    "until_resource_recovery",
]

indicators = {
    "Florida_cultures",
    "Florida_cultures",
    "Florida_cultures",
    "Florida_cultures",
    "Tierra del Fuego (Isla Grande)_cultures",
    "Tierra del Fuego (Isla Grande)_cultures",
    "Haida Nation islands_cultures",
    "Haida Nation islands_cultures",
    "Amazonas_arrival",
    "Haida Nation islands_pop",
    "Alaska_pop",
    "Baja California Sur_pop",
    "California_pop",
    "Amazonas_pop",
}

data["bad"] = (
    (data["Florida_cultures"] < 2) * 1
    + (data["Florida_cultures"] > 20) * 1
    + (data["Florida_cultures"] > 50) * 1
    + (data["Florida_cultures"] > 100) * 1
    + (data["Tierra del Fuego (Isla Grande)_cultures"] < 2) * 1
    + (data["Tierra del Fuego (Isla Grande)_cultures"] > 20) * 1
    + (data["Haida Nation islands_cultures"] > 2) * 1
    + (data["Haida Nation islands_cultures"] > 10) * 1
    + (data["Amazonas_arrival"] < 400) * 1
    + (data["Amazonas_arrival"] < 600) * 1
    + (data["Haida Nation islands_pop"] < 0.25) * 1
    + (data["Alaska_pop"] < 0.1) * 1
    + (data["Baja California Sur_pop"] < 0.3) * 1
    + (data["California_pop"] < 0.6) * 1
    + (data["Amazonas_pop"] < 0.6) * 1
)

x = []
for i, indicator in enumerate(indicators):
    for p, predictor in enumerate(all_predictors):
        a = plt.subplot(
            len(indicators),
            len(all_predictors),
            len(all_predictors) * i + p + 1,
            sharex=x[p] if i else None,
            sharey=plt.gca() if p else None,
        )
        if "_pop" not in indicator:
            plt.yscale("log")
        if i != len(indicators) - 1:
            a.tick_params(labelbottom=False)
        else:
            plt.xlabel(predictor)
        if i == 0:
            x.append(plt.gca())
        if p != 0:
            a.tick_params(labelleft=False)
        else:
            plt.ylabel(indicator)
        plt.scatter(
            data[predictor], data[indicator], c=data["bad"], alpha=0.2, marker="+"
        )
plt.subplots_adjust(hspace=0, wspace=0, left=0.035, right=1, top=1, bottom=0.05)
plt.show()

predictors = all_predictors[:]
models = {}
for indicator in indicators:
    model = sklearn.ensemble.RandomForestRegressor()
    model.fit(
        data[predictors][numpy.isfinite(data[indicator])],
        data[indicator][numpy.isfinite(data[indicator])],
        data["end_time"][numpy.isfinite(data[indicator])],
    )
    data["est_" + indicator] = model.predict(data[predictors])
    print(indicator)
    print(dict(zip(predictors, model.feature_importances_)))
    print()
    models[indicator] = model

    plt.scatter(data[indicator], data["est_" + indicator])
    plt.title(indicator)
    plt.show()

parameter_choices = {}
for key in predictors:
    choices = set()
    for value, badness in data[[key, "bad"]].groupby(key):
        choices.add(value)
    print(key)
    print(choices)
    parameter_choices[key] = sorted(choices)

data.to_csv("summary.csv")

central_points = [
    {
        "culture_mutation_rate": 0.004,
        "culture_dimensionality": 63,
        "cooperation_threshold": 7,
        "minimum_adaptation": 0.5,
        "fight_deadliness": 0.09999999986030161,
        "maximum_resources_one_adult_can_harvest": 10000000000.0,
        "enemy_discount": 0.4,
        "season_length_in_years": 0.16666666666666666,
        "until_resource_recovery": 10.0,
    },
    {
        "culture_mutation_rate": 0.004,
        "culture_dimensionality": 62,
        "cooperation_threshold": 7,
        "minimum_adaptation": 0.5,
        "fight_deadliness": 0.09999999986030161,
        "maximum_resources_one_adult_can_harvest": 10000000000.0,
        "enemy_discount": 0.4,
        "season_length_in_years": 0.166666666667,
        "until_resource_recovery": 3.0000000000029994,
    },
    {
        "culture_mutation_rate": 0.004,
        "culture_dimensionality": 64,
        "cooperation_threshold": 8,
        "minimum_adaptation": 0.5,
        "fight_deadliness": 0.09999999986030161,
        "maximum_resources_one_adult_can_harvest": 10000000000.0,
        "enemy_discount": 0.5,
        "season_length_in_years": 0.166666666667,
        "until_resource_recovery": 59.999999999879996,
    },
    {
        "culture_mutation_rate": 0.004,
        "culture_dimensionality": 64,
        "cooperation_threshold": 13,
        "minimum_adaptation": 0.8,
        "fight_deadliness": 0.049999999813735485,
        "maximum_resources_one_adult_can_harvest": 10000000000.0,
        "enemy_discount": 0.5,
        "season_length_in_years": 0.166666666667,
        "until_resource_recovery": 3.0000000000029994,
    },
    {
        "culture_mutation_rate": 0.004,
        "culture_dimensionality": 64,
        "cooperation_threshold": 6,
        "minimum_adaptation": 0.8,
        "fight_deadliness": 0.09999999986030161,
        "maximum_resources_one_adult_can_harvest": 7.0,
        "enemy_discount": 0.7578582832549999,
        "season_length_in_years": 0.16666666666666666,
        "until_resource_recovery": 5.999999999988001,
    },
    {
        "culture_mutation_rate": 0.006,
        "culture_dimensionality": 48,
        "cooperation_threshold": 10,
        "minimum_adaptation": 0.5,
        "fight_deadliness": 0.049999999813735485,
        "maximum_resources_one_adult_can_harvest": 8.0,
        "enemy_discount": 0.3,
        "season_length_in_years": 0.166666666667,
        "until_resource_recovery": 59.999999999879996,
    },
    {
        "culture_mutation_rate": 0.004,
        "culture_dimensionality": 36,
        "cooperation_threshold": 11,
        "minimum_adaptation": 0.8,
        "fight_deadliness": 0.049999999813735485,
        "maximum_resources_one_adult_can_harvest": 10000000000.0,
        "enemy_discount": 0.5,
        "season_length_in_years": 0.16666666666666666,
        "until_resource_recovery": 5.0,
    },
    {
        "culture_mutation_rate": 0.004,
        "culture_dimensionality": 64,
        "cooperation_threshold": 8,
        "minimum_adaptation": 0.5,
        "fight_deadliness": 0.5,
        "maximum_resources_one_adult_can_harvest": 10000000000.0,
        "enemy_discount": 0.5,
        "season_length_in_years": 0.16666666666666666,
        "until_resource_recovery": 5.0,
    },
    {
        "culture_mutation_rate": 0.004,
        "culture_dimensionality": 62,
        "cooperation_threshold": 5,
        "minimum_adaptation": 0.5,
        "fight_deadliness": 0.6666666665114462,
        "maximum_resources_one_adult_can_harvest": 7.0,
        "enemy_discount": 0.7578582832549999,
        "season_length_in_years": 0.166666666667,
        "until_resource_recovery": 10.0,
    },
    {
        "culture_mutation_rate": 0.0045,
        "culture_dimensionality": 64,
        "cooperation_threshold": 11,
        "minimum_adaptation": 0.5,
        "fight_deadliness": 0.049999999813735485,
        "maximum_resources_one_adult_can_harvest": 10000000000.0,
        "enemy_discount": 0.3,
        "season_length_in_years": 0.16666666666666666,
        "until_resource_recovery": 20.0,
    },
    {
        "culture_mutation_rate": 0.006,
        "culture_dimensionality": 62,
        "cooperation_threshold": 5,
        "minimum_adaptation": 0.5,
        "fight_deadliness": 0.049999999813735485,
        "maximum_resources_one_adult_can_harvest": 6.0,
        "enemy_discount": 0.4,
        "season_length_in_years": 0.166666666667,
        "until_resource_recovery": 10.0,
    },
    {
        "culture_mutation_rate": 0.002,
        "culture_dimensionality": 62,
        "cooperation_threshold": 3,
        "minimum_adaptation": 0.9,
        "fight_deadliness": 0.6666666665114462,
        "maximum_resources_one_adult_can_harvest": 10000000000.0,
        "enemy_discount": 0.5,
        "season_length_in_years": 0.166666666667,
        "until_resource_recovery": 3.0,
    },
    {
        "culture_mutation_rate": 0.008,
        "culture_dimensionality": 64,
        "cooperation_threshold": 6,
        "minimum_adaptation": 0.5,
        "fight_deadliness": 0.049999999813735485,
        "maximum_resources_one_adult_can_harvest": 7.0,
        "enemy_discount": 0.5,
        "season_length_in_years": 0.16666666666666666,
        "until_resource_recovery": 10.0,
    },
    {
        "culture_mutation_rate": 0.006,
        "culture_dimensionality": 63,
        "cooperation_threshold": 5,
        "minimum_adaptation": 0.5,
        "fight_deadliness": 0.049999999813735485,
        "maximum_resources_one_adult_can_harvest": 8.0,
        "enemy_discount": 0.3,
        "season_length_in_years": 0.166666666667,
        "until_resource_recovery": 5.999999999987999,
    },
    {
        "culture_mutation_rate": 0.001,
        "culture_dimensionality": 63,
        "cooperation_threshold": 13,
        "minimum_adaptation": 0.5,
        "fight_deadliness": 0.049999999813735485,
        "maximum_resources_one_adult_can_harvest": 8.0,
        "enemy_discount": 0.3,
        "season_length_in_years": 0.166666666667,
        "until_resource_recovery": 3.333333333333333,
    },
    {
        "culture_mutation_rate": 0.0043,
        "culture_dimensionality": 63,
        "cooperation_threshold": 13,
        "minimum_adaptation": 1.0,
        "fight_deadliness": 0.09999999986030161,
        "maximum_resources_one_adult_can_harvest": 10000000000.0,
        "enemy_discount": 0.5,
        "season_length_in_years": 0.166666666667,
        "until_resource_recovery": 10.0,
    },
]

x = []
for i, indicator in enumerate(indicators):
    for p, predictor in enumerate(all_predictors):
        a = plt.subplot(
            len(indicators),
            len(all_predictors),
            len(all_predictors) * i + p + 1,
            sharex=x[p] if i else None,
            sharey=plt.gca() if p else None,
        )
        if "pop" not in indicator:
            plt.yscale("log")
        if i != len(indicators) - 1:
            a.tick_params(labelbottom=False)
        else:
            plt.xlabel(predictor)
        if i == 0:
            x.append(plt.gca())
        if p != 0:
            a.tick_params(labelleft=False)
        else:
            plt.ylabel(indicator)
        for central_point in central_points:
            vary_parameter = [
                [central_point[q] if q != predictor else choice for q in predictors]
                for choice in parameter_choices[predictor]
            ] + [[central_point[q] for q in predictors]]
            predictions = models[indicator].predict(vary_parameter)
            plt.plot(
                parameter_choices[predictor],
                predictions[:-1],
                alpha=1 / len(central_points),
                c="k",
            )
            plt.scatter(central_point[predictor], predictions[-1], c="k", marker="+")
plt.subplots_adjust(hspace=0, wspace=0, left=0.035, right=1, top=1, bottom=0.05)
plt.show()


min = numpy.infty
all_param = [x for x in tqdm(itertools.product(*parameter_choices.values()))]
numpy.random.shuffle(all_param)
try:
    for param in tqdm(all_param):
        parameter_values = dict(zip(parameter_choices, param))
        characteristics = {k: models[k].predict([param])[0] for k in indicators}

        badness = (
            numpy.hstack(
                (
                    (characteristics["Florida_cultures"] < 2),
                    (characteristics["Florida_cultures"] > 20),
                    (characteristics["Florida_cultures"] > 50),
                    (characteristics["Florida_cultures"] > 100),
                    (characteristics["Tierra del Fuego (Isla Grande)_cultures"] < 2),
                    (characteristics["Tierra del Fuego (Isla Grande)_cultures"] > 20),
                    (characteristics["Haida Nation islands_cultures"] > 2),
                    (characteristics["Haida Nation islands_cultures"] > 10),
                    (characteristics["Amazonas_arrival"] < 400),
                    (characteristics["Amazonas_arrival"] < 600),
                    (characteristics["Haida Nation islands_pop"] < 0.4),
                    (characteristics["Haida Nation islands_pop"] < 0.3),
                    (characteristics["Alaska_pop"] < 0.4),
                    (characteristics["Alaska_pop"] < 0.3),
                    (characteristics["Baja California Sur_pop"] < 0.4),
                    (characteristics["Baja California Sur_pop"] < 0.3),
                    (characteristics["California_pop"] < 0.4),
                    (characteristics["California_pop"] < 0.3),
                    (characteristics["Amazonas_pop"] < 0.4),
                    (characteristics["Amazonas_pop"] < 0.3),
                )
            )
            * 1
        )

        if badness.sum() <= min:
            if badness.sum() < min:
                central_points = []
            print(
                " ".join(
                    [
                        "cargo run",
                        "--bin=simulation",
                        "--",
                        "--culture-mutation-rate {culture_mutation_rate}",
                        "--culture-dimensionality {culture_dimensionality}",
                        "--cooperation-threshold {cooperation_threshold}",
                        "--minimum-adaptation {minimum_adaptation}",
                        "--fight-deadliness {fight_deadliness}",
                        "--enemy-discount {enemy_discount}",
                        "--season-length {season_length_in_years}",
                        "--resource-recovery {resource_recovery}",
                        "--harvest-per-person-per-year {maximum_resources_one_adult_can_harvest}",
                    ]
                ).format(
                    resource_recovery=1 / parameter_values["until_resource_recovery"],
                    **parameter_values
                )
            )
            filter = True
            central_points.append(parameter_values)
            for parameter, value in parameter_values.items():
                filter = filter & (data[parameter] == value)
            if filter.any():
                print(data[filter])
            print(characteristics)
            print(badness)
            min = badness.sum()
except KeyboardInterrupt:
    pass

x = []
for i, indicator in enumerate(indicators):
    for p, predictor in enumerate(all_predictors):
        a = plt.subplot(
            len(indicators),
            len(all_predictors),
            len(all_predictors) * i + p + 1,
            sharex=x[p] if i else None,
            sharey=plt.gca() if p else None,
        )
        if i != len(indicators) - 1:
            a.tick_params(labelbottom=False)
        else:
            plt.xlabel(predictor)
        if i == 0:
            x.append(plt.gca())
        if p != 0:
            a.tick_params(labelleft=False)
        else:
            plt.ylabel(indicator)
        if "pop" not in indicator:
            plt.yscale("log")
        for central_point in central_points:
            vary_parameter = [
                [central_point[q] if q != predictor else choice for q in predictors]
                for choice in parameter_choices[predictor]
            ] + [[central_point[q] for q in predictors]]
            predictions = models[indicator].predict(vary_parameter)
            plt.plot(
                parameter_choices[predictor],
                predictions[:-1],
                alpha=1 / len(central_points),
                c="k",
            )
            plt.scatter(central_point[predictor], predictions[-1], c="k", marker="+")
plt.subplots_adjust(hspace=0, wspace=0, left=0.035, right=1, top=1, bottom=0.05)
plt.show()
