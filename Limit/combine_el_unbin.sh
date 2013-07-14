combine -M Asymptotic --minimizerAlgo Minuit2 --minosAlgo stepping -n elHWWlvjj_unbin -m  600 -d  hwwlvj_ggH600_el_unbin.txt  -v 2 
combine -M Asymptotic --minimizerAlgo Minuit2 --minosAlgo stepping -n elHWWlvjj_unbin -m  700 -d  hwwlvj_ggH700_el_unbin.txt  -v 2  
combine -M Asymptotic --minimizerAlgo Minuit2 --minosAlgo stepping -n elHWWlvjj_unbin -m  800 -d  hwwlvj_ggH800_el_unbin.txt  -v 2  
combine -M Asymptotic --minimizerAlgo Minuit2 --minosAlgo stepping -n elHWWlvjj_unbin -m  900 -d  hwwlvj_ggH900_el_unbin.txt  -v 2 
combine -M Asymptotic --minimizerAlgo Minuit2 --minosAlgo stepping -n elHWWlvjj_unbin -m 1000 -d hwwlvj_ggH1000_el_unbin.txt  -v 2 

hadd -f higgisCombin_el_unbin.root higgsCombineelHWWlvjj_unbin.Asymptotic.mH600.root higgsCombineelHWWlvjj_unbin.Asymptotic.mH700.root higgsCombineelHWWlvjj_unbin.Asymptotic.mH800.root higgsCombineelHWWlvjj_unbin.Asymptotic.mH900.root higgsCombineelHWWlvjj_unbin.Asymptotic.mH1000.root 

root -b ../../DrawLimit.C\(\"el\",\"unbin\"\) -q
