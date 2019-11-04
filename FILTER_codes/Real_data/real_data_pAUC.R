
#########################################################################################################

pAUC_compute = function(p_cur, y, pRange=c(1, 0.95), correct=FALSE)
{
	y = as.factor(y)
	rocobj = roc(y, p_cur)
	pAUC = auc(rocobj, partial.auc=pRange, partial.auc.focus="sp", partial.auc.correct=correct)
	pAUC = as.numeric(pAUC)

	return(pAUC)
}

#########################################################################################################


main_path = "/home/FILTER_codes/Real_data/"
if(!dir.exists(main_path)){
	dir.create(main_path, recursive=TRUE)
}
setwd(main_path)

library(pROC)

# m: the number of folds of cross-validation
# tune_type: the criteria for choosing the lambda in FILTER
m = 5
tune_type = "optimal"
file = paste0("Real_result_K0_Ratio100_B100_kmeans_auc_Folds5_", tune_type, ".RData")
load(file)

# pAUCs_all records all the output
pAUCs_all = NULL
#####################################################
## AUC
#####################################################

K = 2
AUC_fusion_2 = mean(index_fusions$'2'[,2])

K = 4
AUC_fusion_4 = mean(index_fusions$'4'[,2])

K = 6
AUC_fusion_6 = mean(index_fusions$'6'[,2])

K = 8
AUC_fusion_8 = mean(index_fusions$'8'[,2])

AUC_cart = mean(index_carts[,2])
AUC_CB = mean(index_CBs[,2])
AUC_rf = mean(index_rfs[,2])
AUC_logistic_refit = mean(index_glms[,2])
AUC_l1_logistic = mean(index_las[,2])
AUC_logistic = mean(index_glms_only[,2])


AUCs = c(AUC_cart, AUC_rf, AUC_CB, AUC_logistic, AUC_logistic_refit, AUC_l1_logistic, 
		AUC_fusion_2, AUC_fusion_4, AUC_fusion_6, AUC_fusion_8)

pAUCs_all = cbind(pAUCs_all, AUCs)



#####################################################
## 1-0.95, correct=FALSE
#####################################################
pRange = c(1, 0.9)
correct = FALSE

K = 2
pAUC_fusion_2 = mean(sapply(1:m, function(i)pAUC_compute(ps_fusion_out$'2'[i,], ys_out[i,], pRange=pRange, correct=correct)))

K = 4
pAUC_fusion_4 = mean(sapply(1:m, function(i)pAUC_compute(ps_fusion_out$'4'[i,], ys_out[i,], pRange=pRange, correct=correct)))

K = 8
pAUC_fusion_8 = mean(sapply(1:m, function(i)pAUC_compute(ps_fusion_out$'8'[i,], ys_out[i,], pRange=pRange, correct=correct)))

K = 6
pAUC_fusion_6 = mean(sapply(1:m, function(i)pAUC_compute(ps_fusion_out$'6'[i,], ys_out[i,], pRange=pRange, correct=correct)))
pAUC_cart = mean(sapply(1:m, function(i)pAUC_compute(ps_cart_out[i,], ys_out[i,], pRange=pRange, correct=correct)))
pAUC_CB = mean(sapply(1:m, function(i)pAUC_compute(ps_CB_out[i,], ys_out[i,], pRange=pRange, correct=correct)))
pAUC_rf = mean(sapply(1:m, function(i)pAUC_compute(ps_rf_out[i,], ys_out[i,], pRange=pRange, correct=correct)))
pAUC_logistic_refit = mean(sapply(1:m, function(i)pAUC_compute(ps_glm_out[i,], ys_out[i,], pRange=pRange, correct=correct)))
pAUC_l1_logistic = mean(sapply(1:m, function(i)pAUC_compute(ps_la_out[i,], ys_out[i,], pRange=pRange, correct=correct)))
pAUC_logistic = mean(sapply(1:m, function(i)pAUC_compute(ps_glm_only_out[i,], ys_out[i,], pRange=pRange, correct=correct)))

pAUCs = c(pAUC_cart, pAUC_rf, pAUC_CB, pAUC_logistic, pAUC_logistic_refit, pAUC_l1_logistic, 
		pAUC_fusion_2, pAUC_fusion_4, pAUC_fusion_6, pAUC_fusion_8)

pAUCs_all = cbind(pAUCs_all, pAUCs)


#####################################################
## 1-0.95, correct=TRUE
#####################################################
pRange = c(1, 0.9)
correct = TRUE

K = 2
pAUC_fusion_2 = mean(sapply(1:m, function(i)pAUC_compute(ps_fusion_out$'2'[i,], ys_out[i,], pRange=pRange, correct=correct)))

K = 4
pAUC_fusion_4 = mean(sapply(1:m, function(i)pAUC_compute(ps_fusion_out$'4'[i,], ys_out[i,], pRange=pRange, correct=correct)))

K = 8
pAUC_fusion_8 = mean(sapply(1:m, function(i)pAUC_compute(ps_fusion_out$'8'[i,], ys_out[i,], pRange=pRange, correct=correct)))

K = 6
pAUC_fusion_6 = mean(sapply(1:m, function(i)pAUC_compute(ps_fusion_out$'6'[i,], ys_out[i,], pRange=pRange, correct=correct)))
pAUC_cart = mean(sapply(1:m, function(i)pAUC_compute(ps_cart_out[i,], ys_out[i,], pRange=pRange, correct=correct)))
pAUC_CB = mean(sapply(1:m, function(i)pAUC_compute(ps_CB_out[i,], ys_out[i,], pRange=pRange, correct=correct)))
pAUC_rf = mean(sapply(1:m, function(i)pAUC_compute(ps_rf_out[i,], ys_out[i,], pRange=pRange, correct=correct)))
pAUC_logistic_refit = mean(sapply(1:m, function(i)pAUC_compute(ps_glm_out[i,], ys_out[i,], pRange=pRange, correct=correct)))
pAUC_l1_logistic = mean(sapply(1:m, function(i)pAUC_compute(ps_la_out[i,], ys_out[i,], pRange=pRange, correct=correct)))
pAUC_logistic = mean(sapply(1:m, function(i)pAUC_compute(ps_glm_only_out[i,], ys_out[i,], pRange=pRange, correct=correct)))

pAUCs = c(pAUC_cart, pAUC_rf, pAUC_CB, pAUC_logistic, pAUC_logistic_refit, pAUC_l1_logistic, 
		pAUC_fusion_2, pAUC_fusion_4, pAUC_fusion_6, pAUC_fusion_8)

pAUCs_all = cbind(pAUCs_all, pAUCs)

rownames(pAUCs_all) = c("CART", "RF", "CB", "logistic", "logistic_refit", "l1_logistic", 
				"FILTER-2", "FILTER-4", "FILTER-6", "FILTER-8")
colnames(pAUCs_all) = c("AUC", "pAUC", "pAUC(standardized)")
write.csv(pAUCs_all, paste0("pAUCs_", tune_type, ".csv"))

