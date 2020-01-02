import os
import numpy as np
import pandas as pd
import codecs

file_dir = os.path.dirname(os.path.abspath(__file__))

############# MAIN PART ########################
if __name__ == "__main__":

    FPS_paths = []
    for root, dirs, files in os.walk(
        os.path.join(file_dir, "benchmark_results_plots", "FPS")
    ):
        for file in files:
            if file.endswith(".png"):
                FPS_paths.append(os.path.join(root, file))

    METHODS_paths = []
    for root, dirs, files in os.walk(
        os.path.join(file_dir, "benchmark_results_plots", "METHODS")
    ):
        for file in files:
            if file.endswith(".png"):
                METHODS_paths.append(os.path.join(root, file))

    md_lines = []
    md_lines.append("# Benchmark Results\n")
    md_lines.append("## Comparing all fingerprints on a given method\n")
    methods = np.unique(
        [str(x).split("/")[-1].split("_")[0] for x in METHODS_paths]
    )
    for method in methods:
        md_lines.append("### Examining {}\n".format(method))
        md_lines.append("#### AUC\n")
        md_lines.append(
            "![](benchmark_results_plots/METHODS/{}_AUC.png)\n\n".format(
                method
            )
        )  # need to nest

        with open(
            "benchmark_results_plots/METHODS/{}_AUC.html".format(method),
            "r",
        ) as f:
            data = f.readlines()
        md_lines.extend(data)

        # auc_html = pd.read_html(
        #     "benchmark_results_plots/METHODS/{}_AUC.html".format(method),
        # )
        # md_lines.append(str(auc_html)+'\n')

        md_lines.append("#### PRC\n")
        md_lines.append(
            "![](benchmark_results_plots/METHODS/{}_PRC.png)\n\n".format(
                method
            )
        )

        # prc_html = pd.read_html(
        #     "benchmark_results_plots/METHODS/{}_PRC.html".format(method),
        # )
        # md_lines.append(str(prc_html)+'\n')

        with open(
            "benchmark_results_plots/METHODS/{}_PRC.html".format(method),
            "r",
        ) as f:
            data = f.readlines()
        md_lines.extend(data)

    md_lines.append("## Comparing all methods on a given fingerprint\n")
    fps = np.unique([str(x).split("/")[-1].split("_")[0] for x in FPS_paths])
    for fp in fps:
        md_lines.append("### Examining {}\n".format(fp))
        md_lines.append("#### AUC\n")
        md_lines.append(
            "![](benchmark_results_plots/FPS/{}_AUC.png)\n\n".format(fp)
        )
        md_lines.append("#### PRC\n")
        md_lines.append(
            "![](benchmark_results_plots/FPS/{}_PRC.png)\n\n".format(fp)
        )

    with open("benchmark_report.md", "w") as f:
        f.writelines(md_lines)

