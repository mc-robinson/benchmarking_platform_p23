#
# $Id$
#
# module to calculate a fingerprint from SMILES

from rdkit import Chem
from rdkit.Chem import MACCSkeys, AllChem
from rdkit.Avalon import pyAvalonTools as fpAvalon
from rdkit.Chem.AtomPairs import Pairs, Torsions
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.Chem.ChemicalFeatures import BuildFeatureFactory
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import DataStructs

import os
import pandas as pd
import numpy as np

# implemented fingerprints:
# ECFC0 (ecfc0), ECFP0 (ecfp0), MACCS (maccs),
# atom pairs (ap), atom pairs bit vector (apbv), topological torsions (tt)
# hashed atom pairs (hashap), hashed topological torsions (hashtt) --> with 1024 bits
# ECFP4 (ecfp4), ECFP6 (ecfp6), ECFC4 (ecfc4), ECFC6 (ecfc6) --> with 1024 bits
# FCFP4 (fcfp4), FCFP6 (fcfp6), FCFC4 (fcfc4), FCFC6 (fcfc6) --> with 1024 bits
# Avalon (avalon) --> with 1024 bits
# long Avalon (laval) --> with 16384 bits
# long ECFP4 (lecfp4), long ECFP6 (lecfp6), long FCFP4 (lfcfp4), long FCFP6 (lfcfp6) --> with 16384 bits
# RDKit with path length = 5 (rdk5), with path length = 6 (rdk6), with path length = 7 (rdk7)
# 2D pharmacophore (pharm) ?????????????
# CAS bitvector

nbits = 1024
n_cas_bits = 7851
longbits = 16384

# TO TEST CUSTOM FP
# Please see the code below with the example of the cas fingerprint
# The fingerprints for the compounds should be put in the 
# `/../../data/` folder. The code below also needs to modified 
# by both adding a function to creat the fingerprint and 
# a dictionary entry to fp_dict

# put dataframes here, should load on imports #
file_path = os.path.dirname(os.path.abspath(__file__))
data_path = file_path + '/../../data/'

chembl_smiles_df = pd.read_pickle(
    os.path.join(data_path, "chembl_smiles_df.pkl")
)
dud_smiles_df = pd.read_pickle(os.path.join(data_path, "dud_smiles_df.pkl"))
dud_missing_smiles_df = pd.read_pickle(
    os.path.join(data_path, "dud_missing_smiles_df.pkl")
)
muv_smiles_df = pd.read_pickle(os.path.join(data_path, "muv_smiles_df.pkl"))

total_df = pd.concat(
    [chembl_smiles_df, muv_smiles_df, dud_smiles_df, dud_missing_smiles_df]
)
del chembl_smiles_df
del dud_smiles_df
del dud_missing_smiles_df
del muv_smiles_df


def create_cas_fp(mol):
    smi = Chem.MolToSmiles(mol)
    # convert to rdkit fp
    try:
        fp_arr = list(
            total_df.loc[total_df["CANONICAL_SMILES"] == smi, "CAS_fp"]
        )[0]
    except:
        raise ValueError("CAS FP not found for: ", smi)
    bitstring = "".join(fp_arr.astype(str))
    rdkit_fp = DataStructs.cDataStructs.CreateFromBitString(bitstring)
    return rdkit_fp


# dictionary
fpdict = {}
fpdict["ecfp0"] = lambda m: AllChem.GetMorganFingerprintAsBitVect(
    m, 0, nBits=nbits
)
fpdict["ecfp2"] = lambda m: AllChem.GetMorganFingerprintAsBitVect(
    m, 1, nBits=nbits
)
fpdict["ecfp4"] = lambda m: AllChem.GetMorganFingerprintAsBitVect(
    m, 2, nBits=nbits
)
fpdict["ecfp6"] = lambda m: AllChem.GetMorganFingerprintAsBitVect(
    m, 3, nBits=nbits
)
fpdict["ecfp6_cas_length"] = lambda m: AllChem.GetMorganFingerprintAsBitVect(
    m, 3, nBits=n_cas_bits
)
fpdict["ecfc0"] = lambda m: AllChem.GetMorganFingerprint(m, 0)
fpdict["ecfc2"] = lambda m: AllChem.GetMorganFingerprint(m, 1)
fpdict["ecfc4"] = lambda m: AllChem.GetMorganFingerprint(m, 2)
fpdict["ecfc6"] = lambda m: AllChem.GetMorganFingerprint(m, 3)
fpdict["fcfp2"] = lambda m: AllChem.GetMorganFingerprintAsBitVect(
    m, 1, useFeatures=True, nBits=nbits
)
fpdict["fcfp4"] = lambda m: AllChem.GetMorganFingerprintAsBitVect(
    m, 2, useFeatures=True, nBits=nbits
)
fpdict["fcfp6"] = lambda m: AllChem.GetMorganFingerprintAsBitVect(
    m, 3, useFeatures=True, nBits=nbits
)
fpdict["fcfc2"] = lambda m: AllChem.GetMorganFingerprint(
    m, 1, useFeatures=True
)
fpdict["fcfc4"] = lambda m: AllChem.GetMorganFingerprint(
    m, 2, useFeatures=True
)
fpdict["fcfc6"] = lambda m: AllChem.GetMorganFingerprint(
    m, 3, useFeatures=True
)
fpdict["lecfp4"] = lambda m: AllChem.GetMorganFingerprintAsBitVect(
    m, 2, nBits=longbits
)
fpdict["lecfp6"] = lambda m: AllChem.GetMorganFingerprintAsBitVect(
    m, 3, nBits=longbits
)
fpdict["lfcfp4"] = lambda m: AllChem.GetMorganFingerprintAsBitVect(
    m, 2, useFeatures=True, nBits=longbits
)
fpdict["lfcfp6"] = lambda m: AllChem.GetMorganFingerprintAsBitVect(
    m, 3, useFeatures=True, nBits=longbits
)
fpdict["maccs"] = lambda m: MACCSkeys.GenMACCSKeys(m)
fpdict["ap"] = lambda m: Pairs.GetAtomPairFingerprint(m)
fpdict["tt"] = lambda m: Torsions.GetTopologicalTorsionFingerprintAsIntVect(m)
fpdict[
    "hashap"
] = lambda m: rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(
    m, nBits=nbits
)
fpdict[
    "hashap_cas_length"
] = lambda m: rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(
    m, nBits=n_cas_bits
)
fpdict[
    "hashtt"
] = lambda m: rdMolDescriptors.GetHashedTopologicalTorsionFingerprintAsBitVect(
    m, nBits=nbits
)
fpdict[
    "hashtt_cas_length"
] = lambda m: rdMolDescriptors.GetHashedTopologicalTorsionFingerprintAsBitVect(
    m, nBits=n_cas_bits
)
fpdict["avalon"] = lambda m: fpAvalon.GetAvalonFP(m, nbits)
fpdict["avalon_cas_length"] = lambda m: fpAvalon.GetAvalonFP(m, n_cas_bits)
fpdict["laval"] = lambda m: fpAvalon.GetAvalonFP(m, longbits)
fpdict["rdk5"] = lambda m: Chem.RDKFingerprint(
    m, maxPath=5, fpSize=nbits, nBitsPerHash=2
)
fpdict["rdk6"] = lambda m: Chem.RDKFingerprint(
    m, maxPath=6, fpSize=nbits, nBitsPerHash=2
)
fpdict["rdk6_cas_length"] = lambda m: Chem.RDKFingerprint(
    m, maxPath=6, fpSize=n_cas_bits, nBitsPerHash=2
)
fpdict["rdk7"] = lambda m: Chem.RDKFingerprint(
    m, maxPath=7, fpSize=nbits, nBitsPerHash=2
)
fpdict["cas"] = lambda m: create_cas_fp(m)

def CalculateFP(fp_name, smiles):
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        raise ValueError(
            "SMILES cannot be converted to a RDKit molecules:", smiles
        )
    if fp_name not in fpdict:
        raise ValueError(
            "Fingerprint name cannot be found:", fp_name
        )
    return fpdict[fp_name](m)
