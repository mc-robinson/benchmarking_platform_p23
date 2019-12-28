# file added by Matt Robinson, matthew.robinson@postera.ai
# This is the master script to run the benchmark
#
# Using this file avoids navigating the confusing file structure
# All of the the file navigation is done from within this file,
# Thus, this file is not particularly pretty...

from optparse import OptionParser
import os
import sys


# import functions for scoring step
sys.path.insert(0, os.getcwd() + "/scoring/")
import scoring_functions as scor

# prepare command-line option parser
usage = "usage: %prog [options] arg"
parser = OptionParser(usage)

parser.add_option(
    "-d",
    "--dataset",
    dest="dataset",
    type="string",
    help="""
    dataset to use for benchmarking (default: I, options are I, II, or filtered)

    Details of each dataset can be found in the original publications listed here:
    https://github.com/rdkit/benchmarking_platform
    """,
)
parser.add_option(
    "-n",
    "--num",
    dest="num",
    type="int",
    metavar="INT",
    help="number of query mols (default=10, options are 5,10,20)",
)
parser.add_option(
    "-f",
    "--fingerprint",
    type="string",
    dest="fp",
    help="""
    fingerprint to use for benchmarking (default is ecfp6, options are:
    ecfp0,
    ecfp2,
    ecfp4,
    ecfp6,
    ecfc0,
    ecfc2,
    ecfc4,
    ecfc6,
    fcfp2,
    fcfp4,
    fcfp6,
    fcfc2,
    fcfc4,
    fcfc6,
    lecfp4,
    lecfp6,
    lfcfp4,
    lfcfp6,
    maccs,
    ap,
    tt,
    hashap,
    hashtt,
    avalon,
    laval,
    rdk5,
    rdk6,
    rdk7,
    cas,)
    NOTE: the cas fingerprint may only be used if the requisite data
    is provided in the `data` folder.

    Details regarding the calculation of each fingerprint are given in `scoring/fingerprint_lib.py`
    """,
)
parser.add_option(
    "-m",
    "--method",
    dest="method",
    type="string",
    metavar="NAME",
    help="""
         NAME of method (either ML algorithm or similarity metric) to use
         (default: Tanimoto similarity, other options are: NB, RF, LR,
         or any of the following similarity metrics:
         Dice, Cosine, Russel, Kulczynski, McConnaughey, Manhattan, RogotGoldberg
         """,
)
parser.add_option(
    "-o",
    "--outdir",
    dest="outdir",
    metavar="PATH",
    help="output PATH of results directory (default: `./benchmark_results`)",
)

############# MAIN PART ########################
if __name__ == "__main__":

    print("cwd: ", os.getcwd())

    # read in command line options
    (options, args) = parser.parse_args()
    # set values based on arguemnts

    DEFAULT_DATASET = "I"
    if options.dataset:
        if options.dataset == "I":
            import configuration_file_I as conf

            scoring_path = os.path.join(os.getcwd(), "scoring", "data_sets_I")
            validation_path = os.path.join(
                os.getcwd(), "validation", "data_sets_I"
            )
            analysis_path = os.path.join(
                os.getcwd(), "analysis", "data_sets_I"
            )
        elif options.dataset == "II":
            import configuration_file_II as conf

            scoring_path = os.path.join(os.getcwd(), "scoring", "data_sets_II")
            validation_path = os.path.join(
                os.getcwd(), "validation", "data_sets_II"
            )
            analysis_path = os.path.join(
                os.getcwd(), "analysis", "data_sets_II"
            )
        elif options.dataset == "filtered":
            import configuration_file_I_filtered as conf

            scoring_path = os.path.join(os.getcwd(), "scoring", "data_sets_I")
            validation_path = os.path.join(
                os.getcwd(), "validation", "data_sets_I"
            )
            analysis_path = os.path.join(
                os.getcwd(), "analysis", "data_sets_I"
            )
        else:
            raise RuntimeError("the dataset option was not recognized")
    else:
        if DEFAULT_DATASET == "I":
            import configuration_file_I as conf

            scoring_path = os.path.join(os.getcwd(), "scoring", "data_sets_I")
            validation_path = os.path.join(
                os.getcwd(), "validation", "data_sets_I"
            )
            analysis_path = os.path.join(
                os.getcwd(), "analysis", "data_sets_I"
            )
        elif DEFAULT_DATASET == "II":
            import configuration_file_II as conf

            scoring_path = os.path.join(os.getcwd(), "scoring", "data_sets_II")
            validation_path = os.path.join(
                os.getcwd(), "validation", "data_sets_II"
            )
            analysis_path = os.path.join(
                os.getcwd(), "analysis", "data_sets_II"
            )

    DEFAULT_NUM_QUERY_MOLS = 10
    if options.num:
        num_query_mols = options.num
    else:
        num_query_mols = DEFAULT_NUM_QUERY_MOLS
    scor.checkQueryMols(num_query_mols, conf.list_num_query_mols)

    DEFAULT_FP = "ecfp6"
    if options.fp:
        fp = options.fp  # error check occurs in fingerprint_lib
    else:
        fp = DEFAULT_FP

    DEFAULT_METHOD = "Tanimoto"
    # create bools to tell if ML or Similarity
    MODE = "ML"
    if options.method:
        if options.method in ["NB", "RF", "LR"]:
            method = options.method
        elif options.method in [
            "Tanimoto",
            "Dice",
            "Cosine",
            "Russel",
            "Kulczynski",
            "McConnaughey",
            "Manhattan",
            "RogotGoldberg",
        ]:
            method = options.method
            MODE = "SIM"
        else:
            raise RuntimeError("the method option was not recognized")
    else:
        method = "Tanimoto"
        MODE = "SIM"

    DEFAULT_OUTDIR = os.path.join(os.getcwd(), "benchmark_results/")
    if options.outdir:
        outdir = options.outdir
    else:
        outdir = DEFAULT_OUTDIR
    if not os.path.exists(outdir):
        os.system("mkdir {path}".format(path=outdir))
    scor.checkPath(outdir, "outdir")

    ### run the required benchmark scripts ###

    # first score
    os.chdir(scoring_path)
    print("cwd: ", os.getcwd())
    if MODE == "ML":
        # need to change dir
        os.system(
            "python calculate_scored_lists_{ml_method}.py -n {ml_num} -f {ml_fp} -s {ml_sim} -o {ml_outpath}".format(
                ml_num=num_query_mols,
                ml_fp=fp,
                ml_sim="Tanimoto",
                ml_outpath=os.path.join(outdir, MODE, "scoring", method, fp),
            )
        )
    elif MODE == "SIM":
        with open(os.path.join(scoring_path, "fp_file.txt"), "w") as f:
            f.write(fp)
        os.system(
            "python calculate_scored_lists.py -n {sim_num} -f ./fp_file.txt -s {sim_method} -o {sim_outpath}".format(
                sim_num=num_query_mols,
                sim_method=method,
                sim_outpath=os.path.join(outdir, MODE, "scoring", method, fp),
            )
        )
        os.system("rm fp_file.txt")

    # then validate
    # need to find a way to implement PRC
    os.chdir(validation_path)
    print("cwd: ", os.getcwd())
    if not os.path.exists(os.path.join(validation_path, "methods_file.txt")):
        with open(os.path.join(validation_path, "methods_file.txt"), "w") as f:
            f.write("AUC\n")
            f.write("PRC")
    os.system(
        "python calculate_validation_methods.py -m methods_file.txt -i {val_input_path} -o {val_output_path}".format(
            val_input_path=os.path.join(outdir, MODE, "scoring", method, fp),
            val_output_path=os.path.join(
                outdir, MODE, "validation", method, fp
            ),
        )
    )


    # then analayze
    os.chdir(analysis_path)
    os.system(
        "python run_analysis.py -i {anal_input_path} -o {anal_output_path}".format(
            anal_input_path=os.path.join(
                outdir, MODE, "validation", method, fp
            ),
            anal_output_path=os.path.join(
                outdir, MODE, "analysis", method, fp
            ),
        )
    )

    # might need to run something else to get an fp summary
    # then add in some code to automatically generate summary plots at the end...

    os.system(
        "python run_fp_summary.py -i {fp_anal_input_path}".format(
            fp_anal_input_path=os.path.join(
                outdir, MODE, "analysis", method, fp
            )
        )
    )
