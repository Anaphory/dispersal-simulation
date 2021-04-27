import numpy
import bisect
import pandas
import sklearn
import itertools
import sklearn.ensemble
from matplotlib import pyplot as plt
from osm import areas

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
    + (data["Amazonas_arrival"] < 500) * 1
    + (data["Haida Nation islands_pop"] < 0.25) * 1
    + (data["Alaska_pop"] < 0.1) * 1
    + (data["Baja California Sur_pop"] < 0.3) * 1
    + (data["California_pop"] < 0.6) * 1
    + (data["Amazonas_pop"] < 0.6) * 1
)

predictors = all_predictors[:]
models = {}
for indicator in indicators:
    model = sklearn.ensemble.RandomForestRegressor()
    model.fit(
        data[predictors][numpy.isfinite(data[indicator])],
        data[indicator][numpy.isfinite(data[indicator])],
        numpy.log(data["end_time"][numpy.isfinite(data[indicator])] + 1),
    )
    data["est_" + indicator] = model.predict(data[predictors])
    print(indicator)
    print(dict(zip(predictors, model.feature_importances_)))
    print()
    models[indicator] = model

for indicator in indicators:
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
    parameter_choices[key] = list(choices)

data.to_csv("summary.csv")

min = numpy.infty
while True:
    parameter_values = {
        parameter: values[numpy.random.randint(len(values))]
        for (parameter, values) in parameter_choices.items()
    }
    param = [[parameter_values[i] for i in predictors]]
    characteristics = {k: models[k].predict(param)[0] for k in indicators}

    badness = (
        numpy.hstack(
            (
                (models["Florida_cultures"].predict(param) < 2),
                (models["Florida_cultures"].predict(param) > 20),
                (models["Florida_cultures"].predict(param) > 50),
                (models["Florida_cultures"].predict(param) > 100),
                (models["Tierra del Fuego (Isla Grande)_cultures"].predict(param) < 2),
                (models["Tierra del Fuego (Isla Grande)_cultures"].predict(param) > 20),
                (models["Haida Nation islands_cultures"].predict(param) > 2),
                (models["Haida Nation islands_cultures"].predict(param) > 10),
                (models["Amazonas_arrival"].predict(param) < 500),
                (models["Haida Nation islands_pop"].predict(param) < 0.4),
                (models["Haida Nation islands_pop"].predict(param) < 0.3),
                (models["Alaska_pop"].predict(param) < 0.4),
                (models["Alaska_pop"].predict(param) < 0.3),
                (models["Baja California Sur_pop"].predict(param) < 0.4),
                (models["Baja California Sur_pop"].predict(param) < 0.3),
                (models["California_pop"].predict(param) < 0.4),
                (models["California_pop"].predict(param) < 0.3),
                (models["Amazonas_pop"].predict(param) < 0.4),
                (models["Amazonas_pop"].predict(param) < 0.3),
            )
        )
        * 1
    )

    if badness.sum() <= min:
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
                    "--fight-deadliness {deadliness}",
                    "--enemy-discount {enemy_discount}",
                    "--season-length {season_length_in_years}",
                    "--resource-recovery {resource_recovery}",
                    "--harvest-per-person-per-year {maximum_resources_one_adult_can_harvest}",
                ]
            ).format(
                resource_recovery=1 / parameter_values["until_resource_recovery"],
                deadliness=int(parameter_values["fight_deadliness"] * 2 ** 32 + 0.5),
                **parameter_values
            )
        )
        filter = True
        for parameter, value in parameter_values.items():
            filter = filter & (data[parameter] == value)
        if len(data[filter]):
            print(data[filter])
        print(characteristics)
        print(badness)
        min = badness.sum()
