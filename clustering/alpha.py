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
    prog="DBSCAN_ALPHA_ONLY",
    description='Provide your pdb file and get degree of association')
parser.add_argument('-i', '--input', type=Path, required=True)      # option that takes a value
parser.add_argument('-s', '--step',  type=int, required=False, default=1000, help= "Step for collecting the dataframes. Default: 1000 -> quite big number to take only single last dataframe.")
args = parser.parse_args()

INPUT = args.input
OUTPUT_LF = INPUT.parent / f"{INPUT.name[:-4]}_last_frame{INPUT.name[-4:]}" # the same name with modification
STEP_STORE = args.step

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
    array_lengths = np.array(list(c.values()))/21
    # print(array_lengths)
    # print(np.partition(array_lengths, -2))
    # print(np.partition(array_lengths, -2)[-3:])
    longest_olig = np.sum(np.partition(array_lengths, -2)[-4:])

    # print(f"Longest: {longest_olig}")

    alpha = longest_olig / 100

    return alpha


# get the lastest model
rg_out = subprocess.getoutput(f"rg 'MODEL' {INPUT}")
last_model = rg_out.split()[-1:]
# print(last_model)
models_orig = [int(line.split()[-1]) for line in rg_out.splitlines()]
models = models_orig[::-STEP_STORE] # negative -> huge step take only one value

# take all lines
with open(f"{INPUT}") as infile: lines = [line.strip() for line in infile.readlines()]

xyz_ms       = []
xyz_s        = []
write_flag   = False
collect_flag = False
with open(f"{OUTPUT_LF}", "w") as outfile:
    for idx, line in enumerate(lines):
        if idx == 0: continue  # skip the header
        if idx == 1: outfile.write(f"{line}\n") # we always need to put "CRYST1 ..." line

        if "MODEL" in line:
            # find out the number if this is last one
            _, trial_number = line.split()
            if int(trial_number) in models:
                # starting from the next line we need to write it down!
                collect_flag = True
                # continue
            if int(trial_number) == int(last_model[0]):
                # starting from the next line we need to write it down!
                write_flag = True
                # continue
            continue

        # # stop appending
        if "ENDMDL" in line:
            collect_flag = False
            write_flag = False
            # print("====")
            if len(xyz_s) > 0:
                xyz_ms.append(np.array(xyz_s))
                xyz_s.clear()
            continue

        if "TER" in line: continue
        if "CONECT" in line: continue
        if "END" in line: continue

        if write_flag: outfile.write(f"{line}\n")
        if collect_flag:
            # DBG
            # print(line)
            lc = line.split()
            x, y, z = float(lc[6]), float(lc[7]), float(lc[8])
            # print(f"[ {x}, {y}, {z} ]")
            xyz_s.append([x, y, z])

# many
alphas = []
for xyz in xyz_ms:
    X = xyz[:2100]
    db = DBSCAN(eps=4.0, min_samples=1).fit(X)
    labels = db.labels_
    counter = Counter(labels)
    alpha = degree_of_aggregation(counter)
    # print(f"{INPUT}: {INPUT.name.split('_')[0]} {alpha}")
    alphas.append(alpha)
print(f"[{len(alphas)}] Mean: {np.mean(alphas)} ± {np.std(alphas)}")

# single file
# # X Y Z column of position of first 21 * 100 monomers
# # X = np.genfromtxt("./new_1frame_100mon_and_water.pdb", skip_header=1, skip_footer=1, usecols=[6, 7, 8])
# X = np.genfromtxt(OUTPUT_LF, skip_header=1, skip_footer=1, usecols=[6, 7, 8])
# X = X[:2100] # 21* 100 mon
#
# db = DBSCAN(eps=4.0, min_samples=1).fit(X)
# labels = db.labels_
#
# # print(len(labels))
# # print(labels)
# counter = Counter(labels)
#
# # Number of clusters in labels, ignoring noise if present.
# # n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
# # n_noise_    = list(labels).count(-1)
#
# # print("Estimated number of clusters: %d" % n_clusters_)
# # print("Estimated number of noise points: %d" % n_noise_)
#
# alpha = degree_of_aggregation(counter)
# print(f"{INPUT}: {INPUT.name.split('_')[0]} {alpha}")
#
