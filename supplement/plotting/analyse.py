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
        plt.scatter(data[predictor], data[indicator], c=data["bad"], alpha=0.2, marker="+")
plt.subplots_adjust(hspace=0, wspace=0, left=0.035, right=1, top=1, bottom=0.05)
plt.show()

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
all_param = [x for x in tqdm(itertools.product(*parameter_choices.values()))]
numpy.random.shuffle(all_param)
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
        for parameter, value in parameter_values.items():
            filter = filter & (data[parameter] == value)
        if len(data[filter]):
            print(data[filter])
        print(characteristics)
        print(badness)
        min = badness.sum()
