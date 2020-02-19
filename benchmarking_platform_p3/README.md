Benchmarking Platform
=====================
original version presented in 
S. Riniker, G. Landrum, J. Cheminf., 5, 26 (2013),
DOI: 10.1186/1758-2946-5-26,
URL: http://www.jcheminf.com/content/5/1/26

extended version presented in
S. Riniker, N. Fechner, G. Landrum, J. Chem. Inf. Model., 53, 2829, (2013),
DOI: 10.1021/ci400466r,
URL: http://pubs.acs.org/doi/abs/10.1021/ci400466r

GENERAL USAGE NOTES
-------------------
The virtual-screening process implemented by the benchmarking
platform is divided into three steps:

1) Scoring

2) Validation

3) Analysis

The three steps are run separately and read in the output of the
previous step. In the scoring step, the data from the directories
compounds and query_lists is read in.

The directory compounds contains lists of compounds for 118 targets
from three public data sources: MUV, DUD and ChEMBL. The compound
lists contain the external ID, the internal ID and the SMILES of
each compound.

There are three subsets of targets available:

subset I: 
  88 targets from MUV, DUD & ChEMBL described in J. Cheminf., 5, 26 (2013)
  
subset I filtered: 
  69 targets from MUV, DUD & ChEMBL filtered for difficulty
  described in JCIM (2013), online
  
subset II:
  37 targets from ChEMBL designed for a second VS use case
  described in JCIM (2013), online

The directory query_lists contains training lists for each target
with the indices of randomly selected active and inactive molecules.
Training lists with 5, 10 or 20 active molecules are available.
The number of training decoys is 20 % of the decoys for subsets I
and 10 % for subset II.

The scripts are written in Python and use the open-source
cheminformatics library RDKit (www.rdkit.org) and
machine-learning library scikit-learn (www.scikit-learn.org).

Running a script with the option [--help] gives a description of the 
required and optional input parameters of the script.

#### NOTE ABOUT PLOTS ####

If you would like to create the black and gray plots that appear in the manuscript,
please use the following code in `plot_benchmark_results.py`

```
methods = np.unique([c[0] for c in overall_df.columns])
    fps = np.unique([c[-1] for c in overall_df.columns])
    for method in methods:
        for metric in ["AUC", "PRC"]:
            plt.figure(figsize=(16, 12))
            for fp in fps:
                # new option as per Todd's request
                if fp == 'cas':
                    # line_color = 'black'
                    continue
                else:
                    line_color = 'lightgrey'
                plt.errorbar(
                    df["#"],
                    overall_df[(method, metric, fp)],
                    fmt="-o",
                    yerr=overall_df[(method, (metric + "_std"), fp)],
                    label=fp,
                    color=line_color,
                    markerfacecolor=line_color,
                    zorder=-32
                )
            # now plot cas on top 
            plt.errorbar(
                    df["#"],
                    overall_df[(method, metric, 'cas')],
                    fmt="-o",
                    yerr=overall_df[(method, (metric + "_std"), 'cas')],
                    label='cas',
                    color='black',
                    markerfacecolor='black',
                    zorder=3
                )
            plt.xticks(rotation=90)
            plt.legend()
            plt.xlabel("DATASET")
            plt.ylabel(metric + " (w/ SD error bars)")
            plt.title(
                "Comparing all fingerprints with the {} method".format(method)
            )
            plt.savefig(
                os.path.join(
                    plots_dir, "METHODS", "{}_{}.png".format(method, metric)
                )
            )
```
