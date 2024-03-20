# NUMmodel
implementation used for "The unicellular NUM v.0.91: A trait-based plankton model evaluated in two contrasting biogeographic provinces."

The core library is written in Fortran2008 and is interfaced from matlab or R.

### Requirements
The library requires a recent version of matlab (2021 or later). On windows it requires the Matlab MEX module to be installed (Home -> Add-ons -> Get Add-ons -> MATLAB Support for MinGW-w64 C/C++ Compiler); on mac it requires Xcode to be installed. To run global simulation it further requires that the mapping toolbox is installed.  Compiled versions of the library is available for windows (64 bit), linux and osx.  Compiling the library requires a Fortran compiler, e.g., gfortran.  Use the makefile in the Fortran directory. Edit the compiler and flags in the makefile to suit your operating system and compile by writing: `make`.
The water-column setup used for the article requries the use of a transport matrices which must be downloaded separately and placed in the directory TMs. Transport matrices must be downloaded from http://kelvin.earth.ox.ac.uk/spk/Research/TMM/TransportMatrixConfigs (choose MITgcm_2.8deg).

More general info on the NUM model can be found in the README file at the original NUMmodel github: https://github.com/Kenhasteandersen/NUMmodel/releases/tag/v0.91
### Files
The folder includes the core library, and the original matlab and R scripts for simulations. It also includes the files used for this publications.
Due to the number of files genereated in this study they are not all included here. However, they can be reproduced (with a significant computational effort) running the scripts described below. The statistical summeries for all simulations used in this study is saved in the structures saved in * `matlab/processSimulations/processedResults`.
#### files for running the NUMmodel: 
* `matlab/baserunWatercolumn_parameterSensitivityAnalysis.m`.  The script used for running the 200,000 first level simulations for CCE and Station ALOHA
* `matlab/baserunWatercolumn_parameterSensitivityAnalysisReducedParameterSpace.m`. The script used to run the random restricted parameter of the CCE
* `matlab/baserunWatercolumn_parameterSensitivityAnalysisSobolsAnalysis.m`. The script used to run the Sobols Sensitivity analysis on the CCE
These scripts produce a .mat file for each of the simulations. These files are not saved here.
#### Files for processing data
The original result files (.mat) are analyzed and statistics is calculated with the script: 
* `matlab/processSimulations/processResults.m`.  The script processes the result from the NUM model and creates a structure with CRMS, COR, RMSd and STD etc.  
The scripts save the statistics in a structuressaved in * `matlab/processSimulations/processedResults`.
#### Files for plotting figures
* `matlab/processSimulations/PlotScripts.m`.  The script generates most of the figures used in the article. choises has to be made to decide which dataset to plot results from
* `matlab/processSimulations/PlotSobol.m`.  Script that calculates and plots the sobols indexes
* `matlab/processSimulations/TheGood_PlotResult.m`.  Script that plots a subset of results.


