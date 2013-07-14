#!/bin/bash          

FILE='"'$1'"'

echo $FILE
root -l -b -q "TMVAscripts/variables.C($FILE)"
root -l -b -q "TMVAscripts/mvas.C($FILE)"
root -l -b -q "TMVAscripts/mvas.C($FILE,3)"
root -l -b -q "TMVAscripts/efficiencies.C($FILE)"
root -l -b -q "TMVAscripts/correlations.C($FILE)"
#root -l -b -q "TMVAscripts/deviations.C($FILE)"
#root -l -b -q "TMVAscripts/mvasMulticlass.C($FILE)"