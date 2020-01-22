# BENCHMARK CONFIGURATION FILE
File that both describes and specifies the details of the benchamrk run
Simply edit the code-blocks (surrouned with ``` characters) in order to configure your specific parameters of the benchmark.

## Datasets

Options for datasets are datasets I

```
DATASETS = [I]
```

The details of the dataset are below:

#### Dataset I
Dataset I was originally described in 

S. Riniker, G. Landrum, J. Cheminf., 5, 26 (2013),
DOI: 10.1186/1758-2946-5-26,
URL: http://www.jcheminf.com/content/5/1/26

The dataset consists of 88 targets from MUV, DUD, and ChEMBL. 

## Fingerprints
Options for fingerprints are any subset of:
- ecfp0, ecfp2, ecfp4, ecfp6 (1024 bit RDKit Morgan Fingerprint with diameter 0/2/4/6)
- ecfc0, ecfc2, ecfc4, ecfc6 (unhashed RDKit Morgan Fingerprint with diameter 0/2/4/6)
- fcfp2, fcfp4, fcfp6 (1024 bit feature-based RDKit Morgan Fingerprint with diameter 2/4/6)
- fcfc2, fcfc4, fcfc6 (unhashed feature-based RDKit Morgan Fingerprint with diameter 2/4/6)
- lecfp4, lecfp6 (16384 bit RDKit Morgan Fingerprint with diameter 4/6)
- lfcfp4, lfcfp6 (16384 bit feature-based RDKit Morgan Fingerprint with diameter 4/6)
- maccs (MACCS Keys ExplicitBitVect fingerprint)
- ap (Atom Pairs IntSparseIntVect fingerprint)
- hashap (1024 bit hashed Atom Pair fingerprint)
- tt (Topological Torsions LongIntSparseIntVect fingerprint)
- hashtt (1024 bit hashed Topological Torsion fingerprint))
- avalon (1024 bit Avalon Fingerprint)
- laval (16384 bit Avalon Fingerprint)
- rdk5, rdk6, rdk7 (RDKit path-based fingerprint with path-legnth 5/6/7)
- cas (CAS 7851 bit fingerprint. Note that cas data files are needed in `/data` to run cas)
```
FINGERPRINTS = [rdk6, cas]
```

## Methods
Both classical machine learning (ML) methods and basic similarity metrics may be used for activity prediction using the provided fingerprints.

#### ML Methods
The provided ML algorithm options are any subset of 
- Logistic Regression (LR)
- Naive Bayes (NB)
- Random Forests (RF):
```
ML_methods = [LR, NB, RF]
```

The provided similarity metric options are any subset of
- Tanimoto
- Dice
- Cosine
- Russel
- Kulczynski
- McConnaughey
- Manhattan
- RogotGoldberg
```
SIMILARITY_methods = [Tanimoto]
```

## Number of Query Molecules
In each assay, randomly selected subsets of the active and inactive compounds are used for "training" (note that no parameters are fit if using a similarity metric). 

One decision is the number of query active molecules to use for training. Options are 5,10, or 20 active compounds.
```
NUM_QUERY_MOLS = 10
```

