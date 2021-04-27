import pandas
from osm import areas

data = pandas.read_csv(
    "runs_overview.tsv",
    sep="\t",
    names=[
        "file",
        "resource_recovery_per_season",
        "culture_mutation_rate",
        "culture_dimensionality",
        "cooperation_threshold",
        "maximum_resources_one_adult_can_harvest",
        "evidence_needed",
        "payoff_std",
        "minimum_adaptation",
        "fight_deadliness",
        "enemy_discount",
        "season_length_in_years",
        "warfare",
        "end",
        "end_time",
        "mean_pop",
    ]
    + [
        region + "_" + bit for region in areas for bit in ["pop", "cultures", "arrival"]
    ],
)

data["warfare"] = (data["enemy_discount"] != 1.0) | (data["fight_deadliness"] != 0)
data["until_resource_recovery"] = (
    data["season_length_in_years"] / data["resource_recovery_per_season"]
)
data["fight_deadliness"] = data["fight_deadliness"] / 2 ** 32

times = [t for t in data.columns if t.endswith("arrival")]
populations = [t for t in data.columns if t.endswith("pop")]
for t in times:
    data[t] = data[t] * 30 * data["season_length_in_years"]
data["run"] = data["file"].str[-20:]

data = data[data["season_length_in_years"] < 0.2]

data = data.drop_duplicates(subset=["run"], keep='last')

data["finished"] = (data["end"] != "Timeout") * 1 + (data["end"] == "Died out") * 1

data["maximum_resources_one_adult_can_harvest"] = pandas.to_numeric(data["maximum_resources_one_adult_can_harvest"].str.strip().str.split(' ', 1, expand=True)[0])
data["maximum_resources_one_adult_can_harvest"][data["maximum_resources_one_adult_can_harvest"] < 0.5] = 1e50
data["maximum_resources_one_adult_can_harvest"][data["maximum_resources_one_adult_can_harvest"] > 1e10] = 1e10
