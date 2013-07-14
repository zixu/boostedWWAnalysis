combine -M Asymptotic --minimizerAlgo Minuit2 --minosAlgo stepping -n muHWWlvjj_counting -m  600 -d  hwwlvj_ggH600_mu_counting.txt   
combine -M Asymptotic --minimizerAlgo Minuit2 --minosAlgo stepping -n muHWWlvjj_counting -m  700 -d  hwwlvj_ggH700_mu_counting.txt   
combine -M Asymptotic --minimizerAlgo Minuit2 --minosAlgo stepping -n muHWWlvjj_counting -m  800 -d  hwwlvj_ggH800_mu_counting.txt   
combine -M Asymptotic --minimizerAlgo Minuit2 --minosAlgo stepping -n muHWWlvjj_counting -m  900 -d  hwwlvj_ggH900_mu_counting.txt   
combine -M Asymptotic --minimizerAlgo Minuit2 --minosAlgo stepping -n muHWWlvjj_counting -m 1000 -d hwwlvj_ggH1000_mu_counting.txt   

hadd -f higgisCombin_mu_counting.root higgsCombinemuHWWlvjj_counting.Asymptotic.mH600.root  higgsCombinemuHWWlvjj_counting.Asymptotic.mH700.root  higgsCombinemuHWWlvjj_counting.Asymptotic.mH800.root  higgsCombinemuHWWlvjj_counting.Asymptotic.mH900.root  higgsCombinemuHWWlvjj_counting.Asymptotic.mH1000.root

root -b ../../DrawLimit.C\(\"mu\",\"counting\"\) -q
