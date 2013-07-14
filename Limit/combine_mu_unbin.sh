combine -M Asymptotic --minimizerAlgo Minuit2 --minosAlgo stepping -n muHWWlvjj_unbin -m  600 -d  hwwlvj_ggH600_mu_unbin.txt  -v 2 
combine -M Asymptotic --minimizerAlgo Minuit2 --minosAlgo stepping -n muHWWlvjj_unbin -m  700 -d  hwwlvj_ggH700_mu_unbin.txt  -v 2  
combine -M Asymptotic --minimizerAlgo Minuit2 --minosAlgo stepping -n muHWWlvjj_unbin -m  800 -d  hwwlvj_ggH800_mu_unbin.txt  -v 2  
combine -M Asymptotic --minimizerAlgo Minuit2 --minosAlgo stepping -n muHWWlvjj_unbin -m  900 -d  hwwlvj_ggH900_mu_unbin.txt  -v 2 
combine -M Asymptotic --minimizerAlgo Minuit2 --minosAlgo stepping -n muHWWlvjj_unbin -m 1000 -d hwwlvj_ggH1000_mu_unbin.txt  -v 2 

hadd -f higgisCombin_mu_unbin.root higgsCombinemuHWWlvjj_unbin.Asymptotic.mH600.root higgsCombinemuHWWlvjj_unbin.Asymptotic.mH700.root higgsCombinemuHWWlvjj_unbin.Asymptotic.mH800.root higgsCombinemuHWWlvjj_unbin.Asymptotic.mH900.root higgsCombinemuHWWlvjj_unbin.Asymptotic.mH1000.root 

root -b ../../DrawLimit.C\(\"mu\",\"unbin\"\) -q
