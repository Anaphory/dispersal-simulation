import numpy
import pandas
from osm import areas

data = pandas.read_csv(
    "../simulation-results/runs_overview.tsv",
    sep="\t",
    # names=["file", "resource_recovery_per_season", "culture_mutation_rate", "culture_dimensionality", "cooperation_threshold", "maximum_resources_one_adult_can_harvest", "evidence_needed", "payoff_std", "minimum_adaptation", "fight_deadliness", "enemy_discount", "season_length_in_years", "warfare", "end", "end_time", "mean_pop",] + [region + "_" + bit for region in areas for bit in ["pop", "cultures", "arrival"]],
)

data["warfare"] = (data["enemy_discount"] != 1.0) | (data["fight_deadliness"] != 0)
data["until_resource_recovery"] = (
    data["season_length_in_years"] / data["resource_recovery_per_season"]
)
data["fight_deadliness"] = data["fight_deadliness"] / 2 ** 32

data.loc[data["end"] == "Died out", "last"] = 100000

times = [t for t in data.columns if t.endswith("arrival")]
populations = [p for p in data.columns if p.endswith("relative")]
for t in times:
    r = t[:-7] + "relative"
    p = t[:-7] + "persistence"
    # The default logging is every 30 seasons, so rescale some times:
    data[p] = data[p] * 30 * data["season_length_in_years"]
    data[t] = data[t] * 30 * data["season_length_in_years"]
    data.loc[numpy.isnan(data[t]), r] = numpy.nan
    data.loc[data["end"] == "Died out", r] = 0.0
    c = t[:-7] + "cultures"
    data.loc[data["end"] == "Died out", c] = 0
data["run"] = data["Path"].str[-20:]


data = data[data["season_length_in_years"] < 0.2]

data = data.drop_duplicates(subset=["run"], keep='last')

data["finished"] = (data["end"] != "Timeout") * 1 + (data["end"] == "Died out") * 1

data["maximum_resources_one_adult_can_harvest"] = pandas.to_numeric(data["maximum_resources_one_adult_can_harvest"].str.strip().str.split(' ', 1, expand=True)[0])
data.loc[data["maximum_resources_one_adult_can_harvest"] < 0.5, "maximum_resources_one_adult_can_harvest"] = 1e50
data.loc[data["maximum_resources_one_adult_can_harvest"] > 1e10, "maximum_resources_one_adult_can_harvest"] = 1e10

data["enemy_discount"] = data["enemy_discount"].round(7)
data["fight_deadliness"] = data["fight_deadliness"].round(6)
data["until_resource_recovery"] = data["until_resource_recovery"].round(3)
