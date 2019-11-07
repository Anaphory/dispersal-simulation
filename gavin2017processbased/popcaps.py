import csv
import numpy.random

popcaps = []

with open("../binford/data-raw/HG_Output_final.csv") as binford_data:
    for people in csv.DictReader(binford_data, skipinitialspace=True):
        if people["STATE"].startswith("AUSTRALIA"):
            # Exclude AUS populations
            continue
        if int(people["CLIM"]) < 3:
            # Exclude arctic and sub-arctic populations
            continue
        popcaps.append(float(people["TLPOP"]))

def random_popcap():
    return popcaps[numpy.random.randint(len(popcaps))]
