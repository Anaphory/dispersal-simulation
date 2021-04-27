import numpy
import bisect
import operator
import functools

from interpret_data import data, times, populations

all_predictors = [
    "culture_mutation_rate",
    "culture_dimensionality",
    "cooperation_threshold",
    "minimum_adaptation",
    "fight_deadliness",
    "enemy_discount",
    "season_length_in_years",
    "until_resource_recovery",
]

# data = data[data["run"].str.startswith("DA")]

data = data.drop_duplicates(subset=["run"], keep='last')

data["bad"] = (
    (data["end"] == "Died out") * 1 +
    (data["Florida_cultures"] < 2) * 1 +
    (data["Florida_cultures"] > 20) * 1 +
    (data["Florida_cultures"] > 50) * 1 +
    (data["Florida_cultures"] > 100) * 1 +
    (data["Tierra del Fuego (Isla Grande)_cultures"] < 2) * 1 +
    (data["Tierra del Fuego (Isla Grande)_cultures"] > 20) * 1 +
    (data["Haida Nation islands_cultures"] > 2) * 1 +
    (data["Haida Nation islands_cultures"] > 10) * 1 +
    (data["Amazonas_arrival"] < 500) * 1 +
    (data["Haida Nation islands_pop"] < 0.25) * 1 +
    (data["Alaska_pop"] < 0.1) * 1 +
    (data["Baja California Sur_pop"] < 0.3) * 1 +
    (data["California_pop"] < 0.6) * 1 +
    (data["Amazonas_pop"] < 0.6) * 1
)

data.to_csv("summary.csv")

all_predictors = [
    "culture_mutation_rate",
    "culture_dimensionality",
    "cooperation_threshold",
    "minimum_adaptation",
    "fight_deadliness",
    "enemy_discount",
    "season_length_in_years",
    "until_resource_recovery",
]

values_and_their_fit = {}
for key in all_predictors:
    choices = {}
    for value, badness in data[[key, "bad"]].groupby(key):
        weight = 1 / badness["bad"].min() ** 2
        choices[value] = weight
    print(key)
    print(choices)
    values_and_their_fit[key] = (list(choices.keys()), numpy.cumsum(list(choices.values())))


for i in range(100):
    random = numpy.random.random(size=len(all_predictors))
    parameter_values = {
        key: values[bisect.bisect(weights, weights[-1] * r)]
        for (key, (values, weights)), r in zip(values_and_their_fit.items(), random)
    }
    filter = functools.reduce(operator.and_, (data[p] == v for p, v in parameter_values.items()), True)
    print(parameter_values)
    if filter.any():
        print(data[filter])
