# program to plot the results of the benchmark
# Matt Robinson, matthew.robinson@postera.ai

import os
from glob import glob

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from rdkit import Chem
from rdkit.Chem import AllChem

file_path = os.path.dirname(os.path.abspath(__file__))

ML_folder_paths = [p.strip('/') for p in 
                   glob(os.path.join(file_path, 'ML', 'analysis','*/'))]
