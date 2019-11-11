# TactileSearch
Raw data (.mat, .xls) and code (.m, .R) for analysis of behavioral data from tactile search experiment.

Analysis steps are explained in order in the top-level analysis script TactileSearchAnalysis.m

R directory in /Code/ contains the full phase-aligned and phase-scrambled datasets that were exported to .xls files and the associated R code used to fit linear mixed effects models to the 2 datasets.

The 4 .mat files in the main directory (e.g. sAll_10T_PhSc_FixLinOMP_10000iter.mat) are the results from the time-consuming matching pursuit procedure we ran to fit the two-site model. These files are loaded and analyzed by functions in TactileSearchAnalysis.

Some figures employ [subtightplot](https://www.mathworks.com/matlabcentral/fileexchange/39664-subtightplot) by Felipe G. Nievinski. Downloaded version is included in /Code/.
