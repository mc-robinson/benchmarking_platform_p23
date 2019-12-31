# file added by Matt Robinson, matthew.robinson@postera.ai
# This is the master script to run the benchmarks
#
# Using this file avoids navigating the confusing file structure
# All of the the file navigation is done by calling
# benchmarking_platform_p3/run_benchmark.py

import os
import sys

file_dir = os.path.dirname(os.path.abspath(__file__))

############# MAIN PART ########################
if __name__ == "__main__":

    with open(os.path.join(file_dir, "BENCHMARK_CONFIG.md"), "r") as f:
        config_data = f.readlines()

    for idx, line in enumerate(config_data):
        if idx >= (len(config_data) - 1):
            break
        if ("```" in config_data[idx - 1]) and ("```" in config_data[idx + 1]):
            if "DATASETS" in line:
                ds = [
                    x.strip()
                    for x in line.split("[")[-1].split("]")[0].split(",")
                ]
            if "FINGERPRINTS" in line:
                fps = [
                    x.strip()
                    for x in line.split("[")[-1].split("]")[0].split(",")
                ]
            if "ML_methods" in line:
                ML_methods = [
                    x.strip()
                    for x in line.split("[")[-1].split("]")[0].split(",")
                ]
            if "SIMILARITY_methods" in line:
                SIMILARITY_methods = [
                    x.strip()
                    for x in line.split("[")[-1].split("]")[0].split(",")
                ]
            if "NUM_QUERY_MOLS" in line:
                num_query_mols = int(line.split("=")[-1].strip())

    # run in the subdirectory then move to this directory
    os.chdir(os.path.join(file_dir, "benchmarking_platform_p3"))

    for d in ds:
        for fp in fps:
            for ML_method in ML_methods:
                os.system(
                    "python run_benchmark.py -m {ml_method} -n {ml_num} -f {ml_fp} -d {ml_dataset}".format(
                        ml_method=ML_method,
                        ml_num=num_query_mols,
                        ml_fp=fp,
                        ml_dataset=d,
                    )
                )
    
    for d in ds:
        for fp in fps:
            for SIM_method in SIMILARITY_methods:
                os.system(
                    "python run_benchmark.py -m {sim_method} -n {sim_num} -f {sim_fp} -d {sim_dataset}".format(
                        sim_method=SIM_method,
                        sim_num=num_query_mols,
                        sim_fp=fp,
                        sim_dataset=d,
                    )
                )

    os.system("python plot_benchmark_results.py")  
    os.system("mv {} {}".format('./benchmark_results', file_dir))  
    os.system("mv {} {}".format('./benchmark_results_plots', file_dir))       


     

