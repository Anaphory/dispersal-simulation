import numpy
import bisect
import pandas
import sklearn
import argparse
import itertools
import sklearn.ensemble
from matplotlib import pyplot as plt
from osm import areas

from tqdm import tqdm
from interpret_data import data, times, populations

parser = argparse.ArgumentParser()
parser.add_argument("--quick", action="store_true", default=False)
parser.add_argument("--skip", action="store_true", default=False)
args = parser.parse_args()

data = data[
    (0.1666 < data["season_length_in_years"])
    & (data["season_length_in_years"] < 0.1667)
]

all_predictors = [
    "culture_mutation_rate",
    "culture_dimensionality",
    "cooperation_threshold",
    "minimum_adaptation",
    "fight_deadliness",
    "maximum_resources_one_adult_can_harvest",
    "enemy_discount",
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
    "Alaska_relative",
    "Baja California Sur_relative",
    "stdev_logrelpop",
    "Louisiana_persistence",
    "Cuba_persistence",
}


def stdev_logrelpop(data):
    return numpy.ma.masked_invalid(numpy.log(data[populations])).std(axis=1)


data["stdev_logrelpop"] = stdev_logrelpop(data)


def good(data):
    return (
        (data["Florida_cultures"] >= 2) * 1,
        (data["Florida_cultures"] <= 20) * 1,
        (data["Florida_cultures"] <= 50) * 1,
        (data["Florida_cultures"] <= 100) * 1,
        (data["Tierra del Fuego (Isla Grande)_cultures"] >= 2) * 1,
        (data["Tierra del Fuego (Isla Grande)_cultures"] <= 20) * 1,
        (data["Haida Nation islands_cultures"] <= 2) * 1,
        (data["Haida Nation islands_cultures"] <= 10) * 1,
        (data["stdev_logrelpop"] <= 1.0) * 1,
        (data["stdev_logrelpop"] <= 0.9) * 1,
        (data["Alaska_relative"] >= 0.1) * 1,
        (data["Baja California Sur_relative"] >= 0.1) * 1,
        (data["Louisiana_persistence"] >= 120) * 1,
        (data["Cuba_persistence"] >= 120) * 1,
    )


data["good"] = sum(good(data))

# Check which parameter has ben used with which values
parameter_choices = {}
for key in all_predictors:
    choices = set()
    for value, goodness in data[[key, "good"]].groupby(key):
        choices.add(value if int(value) != value else int(value))
    print(key)
    print(choices)
    parameter_choices[key] = sorted(choices)

# Define indicator estimators, with visual control using plots
models = {}
for indicator in indicators:
    model = sklearn.ensemble.RandomForestRegressor()
    model.fit(
        data[all_predictors][numpy.isfinite(data[indicator])],
        data[indicator][numpy.isfinite(data[indicator])],
        data["last"][numpy.isfinite(data[indicator])],
    )
    data["est_" + indicator] = model.predict(data[all_predictors])
    print(indicator)
    print(dict(zip(all_predictors, model.feature_importances_)))
    print()
    models[indicator] = model

    if not args.quick:
        plt.scatter(data[indicator], data["est_" + indicator])
        plt.title(indicator)
        plt.show()

# Sort and store those data points and estimates
data.sort_values(
    ["good", "Baja California Sur_relative"],
    ascending=[False, False],
    na_position="last",
    inplace=True,
)

data.to_csv("summary.csv")

# Plot actual analyses:
# First, assemble all analyses that deviate in at most two parameters from each other.
neighbor_points = {}
for r, row in data.iterrows():
    values = tuple(row[k] for k in all_predictors)
    neighbor_points.setdefault(values, [[], [], []])
for r, row in data.iterrows():
    values = tuple(row[k] for k in all_predictors)
    for comparison in neighbor_points:
        s = (numpy.array(comparison) != numpy.array(values)).sum()
        if s < len(neighbor_points[comparison]):
            neighbor_points[comparison][s].append(r)

best_with_most_data = sorted(
    neighbor_points,
    key=lambda t: (
        data.loc[neighbor_points[t][0]]["good"].max(),
        len(neighbor_points[t][0])
        + 2 * len(neighbor_points[t][1])
        + len(neighbor_points[t][2]),
    ),
    reverse=True,
)

central_points = []
for central_point in best_with_most_data:
    copies = data.loc[neighbor_points[central_point][0]]

    print(central_point, list(copies["good"]), neighbor_points[central_point])

    if copies["good"].max() <= data["good"].max() - 2:
        break
    central_points.append(dict(zip(all_predictors, central_point)))

    if not args.skip:
        deviations = data.loc[neighbor_points[central_point][1]]
        important_points = data.loc[
            neighbor_points[central_point][0] + neighbor_points[central_point][1]
        ]
        supplement = data.loc[neighbor_points[central_point][2]]

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
                if i == 0:
                    x.append(plt.gca())
                if "_relative" not in indicator:
                    plt.yscale("log")
                    miny, maxy = (
                        important_points[indicator].min(),
                        important_points[indicator].max(),
                    )
                    if not numpy.isfinite(miny) or not numpy.isfinite(maxy):
                        continue
                    if miny < 1:
                        miny = 1
                    scale = (maxy / miny) ** 0.05
                    if scale < 1:
                        scale = 1
                    plt.ylim(miny / scale, maxy * scale)
                else:
                    plt.ylim(0, 1)
                minx, maxx = (
                    important_points[predictor].min(),
                    important_points[predictor].max(),
                )
                if not numpy.isfinite(minx) or not numpy.isfinite(maxx):
                    continue
                margin = (maxx - minx) / 20
                plt.xlim(minx - margin, maxx + margin)
                if i != len(indicators) - 1:
                    a.tick_params(labelbottom=False)
                else:
                    plt.xlabel(predictor)
                if p != 0:
                    a.tick_params(labelleft=False)
                else:
                    plt.ylabel(indicator.replace("_", "\n"))
                this_line = deviations[deviations[predictor] != central_point[p]]
                plt.scatter(
                    this_line[predictor],
                    this_line[indicator],
                    c="b",
                    marker="+",
                    alpha=0.3,
                )
                off_center = deviations[deviations[predictor] == central_point[p]]
                plt.scatter(
                    off_center[predictor],
                    off_center[indicator],
                    c="k",
                    marker="+",
                    alpha=0.03,
                )
                off_side = supplement[supplement[predictor] != central_point[p]]
                plt.scatter(
                    off_side[predictor],
                    off_side[indicator],
                    c="k",
                    marker="+",
                    alpha=0.03,
                )
                plt.scatter(copies[predictor], copies[indicator], c="r", marker="+")

        plt.subplots_adjust(hspace=0, wspace=0, left=0.035, right=1, top=1, bottom=0.05)
        plt.show()


def program_call(parameter_values, binary="cargo run --release --bin=simulation --"):
    return " ".join(
        [
            binary,
            "--culture-mutation-rate {culture_mutation_rate}",
            "--culture-dimensionality {culture_dimensionality:d}",
            "--cooperation-threshold {cooperation_threshold:d}",
            "--minimum-adaptation {minimum_adaptation}",
            "--fight-deadliness {fight_deadliness}",
            "--enemy-discount {enemy_discount}",
            "--season-length {season_length_in_years}",
            "--resource-recovery {resource_recovery}",
            "--harvest-per-person-per-year {maximum_resources_one_adult_can_harvest}",
        ]
    ).format(
        resource_recovery=1 / parameter_values["until_resource_recovery"],
        season_length_in_years=1 / 6.0,
        **parameter_values,
    )


if not args.quick:
    central_points_per = [[] for _ in range(14)] + [central_points]
    mins = [0 for _ in range(14)]

    all_param = [x for x in tqdm(itertools.product(*parameter_choices.values()))]
    numpy.random.shuffle(all_param)
    try:
        for param in tqdm(all_param):
            parameter_values = dict(zip(parameter_choices, param))
            characteristics = {k: models[k].predict([param])[0] for k in indicators}

            goodness = good(characteristics)

            echo = False
            for i, (m, g) in enumerate(zip(mins, goodness)):
                if not g:
                    continue
                if sum(goodness) >= mins[i]:
                    echo = True
                    if sum(goodness) > mins[i]:
                        mins[i] = sum(goodness)
                        central_points_per[i] = []
                central_points_per[i].append(parameter_values)
            if echo:
                print(program_call(parameter_values))
                print(characteristics)
                print(goodness)
    except KeyboardInterrupt:
        pass

    x = []
    central_points = {
        tuple(point.values())
        for points in central_points_per
        for point in points
    }
    central_points =[ dict(zip(all_predictors, cond)) for cond in central_points]
    print(len(central_points))
    print()
    for i, indicator in tqdm(enumerate(indicators)):
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
            if "relative" not in indicator:
                plt.yscale("log")
            for central_point in central_points:
                vary_parameter = [
                    [
                        central_point[q] if q != predictor else choice
                        for q in all_predictors
                    ]
                    for choice in parameter_choices[predictor]
                ] + [[central_point[q] for q in all_predictors]]
                predictions = models[indicator].predict(vary_parameter)
                plt.plot(
                    parameter_choices[predictor],
                    predictions[:-1],
                    alpha=1 / len(central_points),
                    c="k",
                )
                plt.scatter(
                    central_point[predictor], predictions[-1], c="k", marker="+"
                )
    plt.subplots_adjust(hspace=0, wspace=0, left=0.035, right=1, top=1, bottom=0.05)
    plt.show()

central_points = sorted(
    central_points,
    key=lambda c: (
        c["culture_dimensionality"] >= 60,
        models["Baja California Sur_relative"].predict(
            [[c[k] for k in all_predictors]]
        ),
    ),
)
i = 0

cpalready = []
for cond in central_points:
    from_back = -1
    central_point = cond[from_back]
    while central_point in cpalready:
        from_back -= 1
        central_point = cond[from_back]
    cpalready.append(central_point)

    k = chr(ord("A") + numpy.random.randint(26)) + chr(
        ord("A") + numpy.random.randint(26)
    )
    run_this = open(f"run-P{k}.sh", "w")
    resume_this = open(f"resume-P{k}.sh", "w")
    run = program_call(central_point, binary="~/data/simulation/simulation")
    print(
        f"sbatch --time=8:00:00 --ntasks=1 --cpus-per-task=8 --output=P{k}{i:03d}.log --partition=generic --wrap '{run} --statefile P{k}{i:03d}.state' --parsable > P{k}{i:03d}.jobid",
        file=run_this,
    )
    print(
        f"sbatch --dependency afterok:`tail -n1 P{k}{i:03d}.jobid` --time=8:00:00 --ntasks=1 --cpus-per-task=8 --output=P{k}r{i:03d}.log --partition=generic --wrap '~/data/simulation/resume --resume-from P{k}{i:03d}.state --statefile P{k}{i:03d}.state >> P{k}{i:03d}.log' --parsable >> P{k}{i:03d}.jobid",
        file=resume_this,
    )

    i += 1
    print(
        f"sbatch --time=8:00:00 --ntasks=1 --cpus-per-task=8 --output=P{k}{i:03d}.log --partition=generic --wrap '{run} --statefile P{k}{i:03d}.state' --parsable > P{k}{i:03d}.jobid",
        file=run_this,
    )
    print(
        f"sbatch --dependency afterok:`tail -n1 P{k}{i:03d}.jobid` --time=8:00:00 --ntasks=1 --cpus-per-task=8 --output=P{k}r{i:03d}.log --partition=generic --wrap '~/data/simulation/resume --resume-from P{k}{i:03d}.state --statefile P{k}{i:03d}.state >> P{k}{i:03d}.log' --parsable >> P{k}{i:03d}.jobid",
        file=resume_this,
    )

    i += 1
    nowar = central_point.copy()
    nowar["enemy_discount"] = 1.0
    nowar["fight_deadliness"] = 1.0
    run = program_call(nowar, binary="~/data/simulation/simulation")
    print(
        f"sbatch --time=8:00:00 --ntasks=1 --cpus-per-task=8 --output=P{k}{i:03d}.log --partition=generic --wrap '{run} --statefile P{k}{i:03d}.state' --parsable > P{k}{i:03d}.jobid",
        file=run_this,
    )
    print(
        f"sbatch --dependency afterok:`tail -n1 P{k}{i:03d}.jobid` --time=8:00:00 --ntasks=1 --cpus-per-task=8 --output=P{k}r{i:03d}.log --partition=generic --wrap '~/data/simulation/resume --resume-from P{k}{i:03d}.state --statefile P{k}{i:03d}.state >> P{k}{i:03d}.log' --parsable >> P{k}{i:03d}.jobid",
        file=resume_this,
    )

    for predictor in all_predictors:
        values = central_point.copy()
        for choice in parameter_choices[predictor]:
            if choice == central_point[predictor]:
                continue
            i += 1
            values[predictor] = choice
            run = program_call(values, binary="~/data/simulation/simulation")
            print(
                f"sbatch --time=8:00:00 --ntasks=1 --cpus-per-task=8 --output=P{k}{i:03d}.log --partition=generic --wrap '{run} --statefile P{k}{i:03d}.state' --parsable > P{k}{i:03d}.jobid",
                file=run_this,
            )
            print(
                f"sbatch --dependency afterany:`tail -n1 P{k}{i:03d}.jobid` --time=8:00:00 --ntasks=1 --cpus-per-task=8 --output=P{k}r{i:03d}.log --partition=generic --wrap '~/data/simulation/resume --resume-from P{k}{i:03d}.state --statefile P{k}{i:03d}.state >> P{k}{i:03d}.log' --parsable >> P{k}{i:03d}.jobid",
                file=resume_this,
            )
