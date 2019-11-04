*fused_coef.R can produce the simulations for the coefficients and cut points esimation and the variable selection performance using FILTER.
*fused_piecewise.R can produce the simulations for prediction performance of different methods for the Model (I) in 4.1.2.
*fused_threshold.R can produce the simulations for prediction performance of different methods for the Model (II) in 4.1.2.
*stat_S1.R can combine the results from fused_coef.R and computing all the indexes as well as creating plots.
*stat_S2.R can combine the results from fused_piecewise.R and fused_threshold.R.
*Directories 'coefs' and 'predictions' contain the outputs of the stat_S1.R and stat_S2.R, which include Table 1, Figure 2, Table 2 and Table3 in the paper. 
*Also, the Monte Carlo experiments in the Supplementary C.1 can be conducted by the above codes.