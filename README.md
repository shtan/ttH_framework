Running instructions:

mkdir bin

First compile ttbar_reco_C package:
cd ttbar_reco_C
make

Then compile main package:
cd ..
make

Before running:
.setup.sh

To run:
./bin/test2 0 10
This runs event 0 through event 10.

If changes are made to ttbar_reco_C package, compile that package before compiling main package again.
