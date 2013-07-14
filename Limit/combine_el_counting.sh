combine -M Asymptotic --minimizerAlgo Minuit2 --minosAlgo stepping -n elHWWlvjj_counting -m  600 -d  hwwlvj_ggH600_el_counting.txt   
combine -M Asymptotic --minimizerAlgo Minuit2 --minosAlgo stepping -n elHWWlvjj_counting -m  700 -d  hwwlvj_ggH700_el_counting.txt   
combine -M Asymptotic --minimizerAlgo Minuit2 --minosAlgo stepping -n elHWWlvjj_counting -m  800 -d  hwwlvj_ggH800_el_counting.txt   
combine -M Asymptotic --minimizerAlgo Minuit2 --minosAlgo stepping -n elHWWlvjj_counting -m  900 -d  hwwlvj_ggH900_el_counting.txt   
combine -M Asymptotic --minimizerAlgo Minuit2 --minosAlgo stepping -n elHWWlvjj_counting -m 1000 -d hwwlvj_ggH1000_el_counting.txt   

hadd -f higgisCombin_el_counting.root higgsCombineelHWWlvjj_counting.Asymptotic.mH600.root  higgsCombineelHWWlvjj_counting.Asymptotic.mH700.root  higgsCombineelHWWlvjj_counting.Asymptotic.mH800.root  higgsCombineelHWWlvjj_counting.Asymptotic.mH900.root  higgsCombineelHWWlvjj_counting.Asymptotic.mH1000.root

root -b ../../DrawLimit.C\(\"el\",\"counting\"\) -q
