# Fingerprint Benchmarking Algorithm
- Matt Robinson, matthew.robinson@postera.ai

## Background ##

This benchmarking platform was originally created by Riniker and Landrum and described in:

S. Riniker, G. Landrum, J. Cheminf., 5, 26 (2013),
DOI: 10.1186/1758-2946-5-26,
URL: http://www.jcheminf.com/content/5/1/26

The dataset consists of 88 targets from MUV, DUD, and ChEMBL. 

The original code can be found at https://github.com/rdkit/benchmarking_platform/blob/master/analysis/analysis_functions.py

## Updating the platform ##

The original platform was written in python 2, and was found to be somehwhat difficult to run. Consequently, we have updated the code to be both python 2 and 3 compatible, for ease of use. Futherrmore, we have changed the structure of the code to facilitate benchmarking without the need for navigating the complicated file structure. 

The original code is contained in he subdirectory `benchmarking_platform_p3`. Most of the code is unchanged excepting changes to allow for python 3 usage. However, several files have been changed more significantly in order to simplify usage. For example, most scripts have been changed to accept abolute -- instead of relative paths -- as parameters.

## Running the code

The easiest way to get started is to clone this GitHub repository. 

#### Dependencies
All code should be run in the provided [conda environment](https://docs.conda.io/projects/conda/en/latest/index.html)

This environment can be installed using 
`conda env create -f environment.yml`

#### Parameters
To modify the specific parameters of the benchmark, please visit the `BENCHMARK_CONFIG.md` file and change the desired parameters.

After this configuration file is appropriately modified, run the benchmark using the `run_benchmark.py` script. Finally, an automatically generated markdown report of the results can be produced using produce report.py.

An overview of these steps is below:
0. Make sure the cas-env conda environment is active
1. Modify BENCHMARK_CONFIG.md as desired
2. `python run benchmark.py`
3. `python produce_report.py`

The generated plots are found in `benchmark_results_plots` and the report is located at `benchmark_report.md`.

#### Testing your own fingerprint

Please see and modify `benchmarking_platform_p3/scoring/fingerprint_lib.py`. The file contains the necessary instructions.

## License

The original benchmarking code was distributed with the following license.s

Copyright (c) 2013, Novartis Institutes for BioMedical Research Inc.
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met: 

     * Redistributions of source code must retain the above copyright 
       notice, this list of conditions and the following disclaimer.
     * Redistributions in binary form must reproduce the above
       copyright notice, this list of conditions and the following 
       disclaimer in the documentation and/or other materials provided 
       with the distribution.
     * Neither the name of Novartis Institutes for BioMedical Research Inc. 
       nor the names of its contributors may be used to endorse or promote 
       products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

