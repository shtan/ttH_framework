This repository contains a C++ implementation of a method to improve object resolutions in events containing two top quarks.

The core package in ttbar_reco_C does the reconstruction fit.
The rest of the code takes Madgraph-generated events, puts them in the correct format, and feeds it to the reconstructor in ttbar_reco_C.
It then outputs the results in text files, and generates plots to compare the reconstructed events with non-reconstructed events.

This repo is adapted from an earlier repo, https://github.com/shtan/Top-Reco-4.
That repo was in turn based on original code from Gala Kaufman and Luke Winstrom.
In this repo I have improved the optimisation process to increase convergence rate,
as well as restructured the code to make it more efficient and less error-prone.
I have also added additional functionality (results analysis, plotting).

As of early 2017, most development has moved to a private repo, where the reconstruction is applied to MC events.
However, this repo is still used to run tests on Madgraph events.

## Setup instructions:

mkdir bin

First compile ttbar_reco_C package:
cd ttbar_reco_C
make

Then compile main package:
cd ..
make

If changes are made to ttbar_reco_C package, compile that package before compiling main package again.

Before running (only need to do once per session):
. setup.sh

## To run:

./bin/test2 firstevent lastevent input outdir task

"task" can be either fit, fitplot, or plot
input: Either the input file (if running fit or fitplot) or input folder (if running plot)
    No ending slash if it's a directory.
outdir: output directory.  Don't put ending slash.

-----------------

Examples:
To do the fit, and output the results in text files:
./bin/test2 0 10 ttbb_h_bbbbdue.lhe ./output/test fit

This runs the fit on event 0 through event 10.  The output is written into text files,
one for each combination of fitstatus, dataset, particle, and momentum component.
The text files will be placed in ./output/test.
Make sure this directory exists before running the code, or nothing will be written.

-----------------

To do the fit, output the results in text files, and also plot histograms that
compare best-gen with smeared-gen:
./bin/test2 0 10 ttbb_h_bbbbdue.lhe ./output/test fitplot

The output text files will be in ./output/test, and the plots will be placed in 
./output/test/plots.
Make sure these two directories exist before running the code, or nothing will be written.

-----------------

To make plots using existing output text files:
./bin/test2 0 1 inputdir ./output/test plot

In this case, firstevent and lastevent are not used.
This will read in all the events in the text files in inputdir, and make plots.
The plots will be placed in ./output/test/plots.
Make sure this directory exists before running.