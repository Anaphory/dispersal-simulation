import numpy
import pandas
import sklearn
import itertools
import sklearn.linear_model
from matplotlib import pyplot as plt
from osm import areas


data = pandas.read_csv("runs_overview.tsv", sep="\t", names=["file", "resource_recovery_per_season", "culture_mutation_rate", "culture_dimensionality", "cooperation_threshold", "maximum_resources_one_adult_can_harvest", "evidence_needed", "payoff_std", "minimum_adaptation", "fight_deadliness", "enemy_discount", "season_length_in_years", "warfare", "end", "end_time", "mean_pop"] + [region + "_" + bit for region in areas for bit in ["pop", "arrival", "cultures"]]

data["warfare"] = (data["enemy_discount"] != 1.0) | (data["fight_deadliness"] != 0)
data["until_resource_recovery"] = data["season_length_in_years"] / data["resource_recovery_per_season"]
data["fight_deadliness"] = data["fight_deadliness"] / 2**32
times = [t for t in data.columns if t.endswith(arrival)]
for t in times:
    data[t] = data[t] * 30 * data["season_length_in_years"]
data["run"] = data["file"].str[38:43]

data = data[data["season_length_in_years"] < 0.2]

all_predictors = ['culture_mutation_rate', 'culture_dimensionality', 'cooperation_threshold', 'minimum_adaptation', 'fight_deadliness', 'enemy_discount', 'season_length_in_years', 'until_resource_recovery']

for i, j in itertools.combinations(all_predictors, 2):
    try:
        new_predictor = f"{i}×{j}"
        # data[new_predictor] = data[i] * data[j]
    except TypeError:
        pass

for name, (Y, filter) in {
        "Died Out": (data["end"] == "Died out", None),
        "Years at timeout (20h)": (data["end_time"], data["end"] == "Timeout"),
        "California_arrival": (data["California_arrival"], numpy.isfinite(data["California_arrival"])),
}.items():
    predictors = all_predictors[:]
    if set(Y) == {True, False}:
        model = sklearn.linear_model.LogisticRegression()
    else:
        model = sklearn.linear_model.LinearRegression()
    while True:
        if filter is None:
            model.fit(data[predictors], Y)
        else:
            model.fit(data[predictors][filter], Y[filter])
        importance = (data[predictors].std() * model.coef_.flat / Y.std())
        pred = pandas.DataFrame({"importance": importance, "coefficients": model.coef_.flat})
        pred["abs_importance"] = abs(pred["importance"])
        pred.sort_values("abs_importance", inplace=True)
        if pred["abs_importance"][0] > 1.0 or len(pred["abs_importance"]) <= 8:
            break
        predictors = list(pred.index)[1:]
    print(name)
    print(model.intercept_)
    print(pred)
    print()

pandas.plotting.scatter_matrix(
    data[all_predictors + times]
)
plt.show()

plt.scatter(data["Amazonas_arrival"], data["Amazonas_pop"], c=(data["end"] == "Ended"))
for run, x, y in zip(data["run"], data["Amazonas_arrival"], data["Amazonas_pop"]):
    if numpy.isfinite(x) and numpy.isfinite(y):
        plt.annotate(run, (x, y), fontsize=12)
plt.xlabel("Amazonas_arrival")
plt.ylabel("Amazonas_pop")
plt.show()

plt.scatter(data["Amazonas_arrival"], data["Tierra del Fuego (Isla Grande)_arrival"], c=(data["end"] == "Ended"))
plt.xlabel("Amazonas_arrival")
plt.ylabel("TdF_arrival")
plt.show()

plt.scatter(data["Paja California Sur_pop"], data["Tierra del Fuego (Isla Grande)_arrival"], c=(data["end"] == "Ended"))
for run, x, y in zip(data["run"], data["BCalS_pop"], data["TdF_arrival"]):
    if numpy.isfinite(x) and numpy.isfinite(y):
        plt.annotate(run, (x, y), fontsize=12)
plt.xlabel("BCalS_pop")
plt.ylabel("TdF_arrival")
plt.show()
