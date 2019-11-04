*files in the directory 'Threshold Effect CI' can generate Fig 1 in the paper.
*real_data.csv is the physical examination/survey data we used.
*real_data_prediction.R can perform cross-validation analysis for our data using different methods including FILTER.
*Real_result_K0_B100_kmeans_auc_Folds5_optimal.RData is the output of the real_data_prediction.R with 5-folds cross-validation procedure, and the other parameters are decleared in the paper or in the program.
*real_data_pAUC.R can perform some calculation on the output of real_data_prediction.R to produce Table 4 in Sec 5.
*pAUCs_optimal.csv is the Table 4 in the paper produced by real_data_pAUC.R.
*score_system.R can produce the FILTER-based risk score in 5.1 for our data..
*score_35.txt is the output scores of all selected variables from score_system.R, which produce the Table 5.
*scoringRule.R can produce the plots of Fig 3, 4 and 5 based on the Murphy diagrams.
*Directories 'insamples' and 'outsamples' contain the Murphy diagrams of insample- and outsample-behaviors of different methods. 
