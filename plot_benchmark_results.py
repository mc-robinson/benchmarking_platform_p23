# program to plot the results of the benchmark
# Matt Robinson, matthew.robinson@postera.ai

import os
from glob import glob

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from rdkit import Chem
from rdkit.Chem import AllChem

file_dir = os.path.dirname(os.path.abspath(__file__))
plots_dir = os.path.join(file_dir, "benchmark_results_plots")
if os.path.exists(plots_dir):
    raise ValueError("PLEASE FIRST REMOVE {}".format(plots_dir))
if not os.path.exists(plots_dir):
    os.system("mkdir {}".format(plots_dir))
    os.system("mkdir {}".format(os.path.join(plots_dir, "METHODS")))
    os.system("mkdir {}".format(os.path.join(plots_dir, "FPS")))

for d in ["I"]:  # ["I", "II"]:
    ML_file_paths = []
    for root, dirs, files in os.walk(
        os.path.join(file_dir, "benchmark_results", d, "ML")
    ):
        for file in files:
            if "summary" in file and file.endswith(".txt"):
                ML_file_paths.append(os.path.join(root, file))

    SIM_file_paths = []
    for root, dirs, files in os.walk(
        os.path.join(file_dir, "benchmark_results", d, "SIM")
    ):
        for file in files:
            if "summary" in file and file.endswith(".txt"):
                SIM_file_paths.append(os.path.join(root, file))

    overall_dict = {}
    for fn in ML_file_paths + SIM_file_paths:
        df = pd.read_csv(fn, sep=" ")
        df = df.rename(columns={"std": "AUC_std", "std.1": "PRC_std"})

        fn = str(fn).strip(".txt")
        fp_name = fn.split("/")[-3]
        method_name = fn.split("/")[-4]

        for col in list(df.columns)[1:-1]:
            overall_dict[(method_name, col, fp_name)] = list(df[col])

    overall_df = pd.DataFrame(overall_dict)

    methods = np.unique([c[0] for c in overall_df.columns])
    fps = np.unique([c[-1] for c in overall_df.columns])
    for method in methods:
        for metric in ["AUC", "PRC"]:
            plt.figure(figsize=(16, 12))
            for fp in fps:
                if fp in "avalon":
                    continue
                plt.errorbar(
                    df["#"],
                    overall_df[(method, metric, fp)],
                    fmt="-o",
                    yerr=overall_df[(method, (metric + "_std"), fp)],
                    label=fp,
                )
            plt.xticks(rotation=90)
            plt.legend()
            plt.xlabel("DATASET")
            plt.ylabel(metric + " (w/ SD error bars)")
            plt.title(
                "Comparing all fingerprints with the {} method".format(fp)
            )
            plt.savefig(
                os.path.join(
                    plots_dir, "METHODS", "{}_{}.png".format(method, metric)
                )
            )

    # now do it for each method on each fingerprint
    methods = np.unique([c[0] for c in overall_df.columns])
    fps = np.unique([c[-1] for c in overall_df.columns])
    for fp in fps:
        if fp in "avalon":
            continue
        for metric in ["AUC", "PRC"]:
            plt.figure(figsize=(16, 12))
            for method in methods:
                plt.errorbar(
                    df["#"],
                    overall_df[(method, metric, fp)],
                    fmt="-o",
                    yerr=overall_df[(method, (metric + "_std"), fp)],
                    label=method,
                )
            plt.xticks(rotation=90)
            plt.legend()
            plt.xlabel("DATASET")
            plt.ylabel(metric + " (w/ SD error bars)")
            plt.title("Comparing all methods on the {} fingerprint".format(fp))

            plt.savefig(
                os.path.join(plots_dir, "FPS", "{}_{}.png".format(fp, metric))
            )

