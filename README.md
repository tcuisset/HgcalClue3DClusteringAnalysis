# Code for testing clustering with test beam data

This setup is a	   modified version of the   setup from: \
https://github.com/rovere/Hgcal_testbeam_analysis_2021

The modifications are required to run on the testbeam ntuples created using the following setup: \
https://github.com/sameasy/TestBeamReconstruction/tree/myTESTS

To run: \
make \
sh runclustering.sh 

Brief description of the scripts:\
Runclustering.cc  : Main analyzer code \
TBNtupleAnalyzer.h : Input tree variables \
CLUEAlgo.h : clue2D \
CLUE3DAlgo.h : clue3D (in x-y-L space at this point)\

Produces tree with hits and 2D-3D cluster properties which can be plotted using notebooks inside plotter folder \

Requires uproot, numpy 