#!/usr/bin/env python3

from pathlib import Path
from collections import Counter
import numpy as np
from sklearn.cluster import DBSCAN
import matplotlib.pyplot as plt
from coloring import prepare_VMD_colors
from shutil import rmtree
from shlex import shlex, join
import subprocess

import argparse

parser = argparse.ArgumentParser(
    prog="DBSCAN",
    description='Provide your pdb file and get colored one')
parser.add_argument('-i', '--input', type=Path, required=True)      # option that takes a value
# parser.add_argument('-s', '--step',  type=int, required=False, default=1000, help= "Step for collecting the dataframes. Default: 1000 -> quite big number to take only single last dataframe.")
args = parser.parse_args()

# INPUT = Path("./input.pdb")
INPUT = args.input
OUTPUT_LF = INPUT.parent / f"{INPUT.name[:-4]}_last_frame{INPUT.name[-4:]}" # the same name with modification
# STEP_STORE = args.step

OUTPUT = INPUT.parent / f"{INPUT.name[:-4]}_colored{INPUT.name[-4:]}" # the same name with modification
COLORS_OUTPUT = Path("./colors.idx")
COLORS_RED_OUTPUT = Path("./colors.idx.red")

PLOT          = False
CLEAN_SUCCESS = True


def degree_of_aggregation(c: Counter):
    '''
    Input: Counter(clusters) -- I used it because it is a map:
    cluster ID vs. number of atoms (will get number of monomers if divided by 21)

    Idea:
     *to find the longest polymer and consider it as agregated in the ratio:
                    number of aggregated monomers
           alpha =  ------------------------------
                       total number of monomers

    alpha -- degree of association

    '''
    alpha = 0.0

    # check if the number of monomers are equal to 100 (hardcoded number)
    assert int(sum(c.values())/21) == 100, "Number of monomer is not 100"

    # find the longest oligomer
    longest_olig = max(c.values()) / 21
    print(f"Longest: {longest_olig}")

    alpha = longest_olig / 100

    return alpha


# get the lastest model
rg_out = subprocess.getoutput(f"rg 'MODEL' {INPUT}")
last_model = rg_out.split()[-1:]
# models_orig = [int(line.split()[-1]) for line in rg_out.splitlines()]
# models = models_orig[::-STEP_STORE] # negative -> huge step take only one value
# print(last_model)

# take all lines
with open(f"{INPUT}") as infile:
    lines = [l.strip() for l in infile.readlines()]

write_flag = False
with open(f"{OUTPUT_LF}", "w") as outfile:
    for idx, line in enumerate(lines):
        if idx == 0: continue  # skip the header
        if idx == 1: outfile.write(f"{line}\n") # we always need to put "CRYST1 ..." line

        if "MODEL" in line:
            # find out the number if this is last one
            _, trial_number = line.split()
            # strart appending
            # if int(trial_number) in models:
            if int(trial_number) == int(last_model[0]):
                # starting from the next line we need to write it down!
                write_flag = True
                continue

        # # stop appending
        if "ENDMDL" in line: continue
        # if "ENDMDL" in line:
        #     write_flag = False
        #     continue

        if "TER" in line: continue
        if "CONECT" in line: continue
        if "END" in line: continue

        if write_flag: outfile.write(f"{line}\n")


# X Y Z column of position of first 21 * 100 monomers
# X = np.genfromtxt("./new_1frame_100mon_and_water.pdb", skip_header=1, skip_footer=1, usecols=[6, 7, 8])
X = np.genfromtxt(OUTPUT_LF, skip_header=1, skip_footer=1, usecols=[6, 7, 8])
X = X[:2100] # 21* 100 mon

db = DBSCAN(eps=4.0, min_samples=1).fit(X)
labels = db.labels_

print(len(labels))
print(labels)
counter = Counter(labels)
clusters = []
for key, val in counter.items():
    print(f"[ID] {key} : {val/21}")
    clusters.append(int(val/21))

with open(COLORS_OUTPUT, "w") as colorfile:
    for label in labels: colorfile.write(f"{label}\n")

# HACK FOR VMD Coloring (VMD has only 32 default colors) so I am readucing the number of labels
# in paper we can mention that each color corresponds to new domain
new_colors = []
for label in labels:
    if label > 32: new_label = label - 32
    else:          new_label = label
    new_colors.append(new_label)
with open(COLORS_RED_OUTPUT, "w") as colorfile:
    for label in new_colors: colorfile.write(f"{label}\n")
# END HACK

# Number of clusters in labels, ignoring noise if present.
n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
n_noise_    = list(labels).count(-1)

print("Estimated number of clusters: %d" % n_clusters_)
print("Estimated number of noise points: %d" % n_noise_)

# take the clusters and calculate the distribution on size of clusters
possible_cluster_lengths = np.arange(100) # 0, 1, 2, 3...n
hist = np.histogram(np.array(clusters), bins=possible_cluster_lengths)
print(hist)
assert sum(hist[0]) == n_clusters_, f"{sum(hist[0])} != {n_clusters_}"
print("Check passed!")

# Coloring for VMD
with open(OUTPUT_LF) as infile: init_lines = [il.strip() for il in infile.readlines()]
prepare_VMD_colors(init_lines, colors_input=COLORS_OUTPUT, output=OUTPUT)


# alpha = degree_of_aggregation(counter)
# print(f"Degree of association: {alpha}")

try:
    p = subprocess.Popen(f"vmd_draw {OUTPUT}", shell= True).wait()
except Exception as e:
    print(f"Problem: {e}")

if PLOT:

    unique_labels = set(labels)
    core_samples_mask = np.zeros_like(labels, dtype=bool)
    core_samples_mask[db.core_sample_indices_] = True

    fig = plt.figure(figsize=(12, 12))
    ax = fig.add_subplot(111, projection="3d")
    colors = [plt.cm.Spectral(each) for each in np.linspace(0, 1, len(unique_labels))]
    for k, col in zip(unique_labels, colors):
        if k == -1:
            # Black used for noise.
            col = [0, 0, 0, 1]

        class_member_mask = labels == k

        xy = X[class_member_mask & core_samples_mask]
        ax.plot(
            xy[:, 0],
            xy[:, 1],
            xy[:, 2],
            "o",
            markerfacecolor=tuple(col),
            markeredgecolor="k",
            markersize=14,
        )

        xy = X[class_member_mask & ~core_samples_mask]
        ax.plot(
            xy[:, 0],
            xy[:, 1],
            xy[:, 2],
            "o",
            markerfacecolor=tuple(col),
            markeredgecolor="k",
            markersize=6,
        )

    plt.title(f"Estimated number of clusters: {n_clusters_}")
    plt.show()

if CLEAN_SUCCESS:
    # remove info: colors and colors reduced
    COLORS_OUTPUT.unlink()
    COLORS_RED_OUTPUT.unlink()


