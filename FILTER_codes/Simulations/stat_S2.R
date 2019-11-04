

####################################################################################################

## this version is the negative Brier score in 2007 JASA paper
Brier = function(y, p){
	return(1 + p^2 + (1-p)^2 - 2*dbinom(y, 1, p))
}

compute_scores = function(y, p)
{
	scores_logs = -logs_binom(y=y, size=1, prob=p)
	scores_crps = -crps_binom(y=y, size=1, prob=p)
	scores_brier = -Brier(y=y, p=p)

	results = cbind(scores_logs, scores_crps, scores_brier)
	colnames(results) = c("Logs", "CRPS", "Brier")
	Inf_index = which(scores_logs==-Inf | scores_crps==-Inf | scores_brier==-Inf)
	if(length(Inf_index)>0)	results = results[-Inf_index,]
	
	return(results)
}

####################################################################################################


##############################################################################
## output the indexes and plots in 4.1.2
##############################################################################


library(scoringRules)

# set the working directory
main_path = "/home/FILTER_codes/Simulations/predictions/"
if(!dir.exists(main_path)){
	dir.create(main_path, recursive=TRUE)
}
setwd(main_path)

# parameters should match the settings in the fused_threshold.R and fused_piecewise.R
ns = c(400, 800)
rhos = c(0, 0.5)
p = 500
p0s = c(5, 10)
########################################################################################
# model = "thresh"  is for Model (I) in 4.1.2
# change this parameter to 'piece' if you want to output the Model (II) in 4.1.2
########################################################################################
model = "piece" 
cur_path = paste0(main_path, model, "_results/")
merge = FALSE

out_matrix_all = NULL

for(p0 in p0s)
{
	for(rho in rhos)
	{	
		for(n in ns)
		{
			# set current directory
			setwd(cur_path)
			
			# merge the results if needed
			if(merge){
				save_list = c("index_carts","index_CBs","index_fusions","index_glms",
						"index_glms_only","index_las","index_rfs","ps_cart",
						"ps_CB","ps_fusion","ps_glm","ps_glm_only","ps_la","ps_rf","ys")

				index_fusions_all = NULL
				index_carts_all = NULL
				index_CBs_all = NULL
				index_glms_all = NULL
				index_glms_only_all = NULL
				index_las_all = NULL
				index_rfs_all = NULL
				
				ps_fusion_all = NULL
				ps_cart_all = NULL
				ps_CB_all = NULL
				ps_glm_all = NULL
				ps_glm_only_all = NULL
				ps_la_all = NULL
				ps_rf_all = NULL

				ys_all = NULL
			
				file_format = paste0("n", n, "_rho", rho*10, "_p0", p0, "_p", p, "_loop")
				files = dir(pattern=paste0(file_format, "[0-9]+.RData"))
				for(rfile in files)
				{
					print(rfile)
					load(rfile)
					index_fusions_all = rbind(index_fusions_all, index_fusions)
					index_carts_all = rbind(index_carts_all, index_carts)
					index_CBs_all = rbind(index_CBs_all, index_CBs)
					index_glms_all = rbind(index_glms_all, index_glms)
					index_glms_only_all = rbind(index_glms_only_all, index_glms_only)
					index_las_all = rbind(index_las_all, index_las)
					index_rfs_all = rbind(index_rfs_all, index_rfs)
				
					ps_fusion_all = rbind(ps_fusion_all, ps_fusion)
					ps_cart_all = rbind(ps_cart_all, ps_cart)
					ps_CB_all = rbind(ps_CB_all, ps_CB)
					ps_glm_all = rbind(ps_glm_all, ps_glm)
					ps_glm_only_all = rbind(ps_glm_only_all, ps_glm_only)
					ps_la_all = rbind(ps_la_all, ps_la)
					ps_rf_all = rbind(ps_rf_all, ps_rf)
				
					ys_all = rbind(ys_all, ys)
					file.remove(rfile)
				}

				index_fusions = index_fusions_all
				index_carts = index_carts_all
				index_CBs = index_CBs_all
				index_glms = index_glms_all
				index_glms_only = index_glms_only_all
				index_las = index_las_all
				index_rfs = index_rfs_all
				
				ps_fusion = ps_fusion_all
				ps_cart = ps_cart_all
				ps_CB = ps_CB_all
				ps_glm = ps_glm_all
				ps_glm_only = ps_glm_only_all
				ps_la = ps_la_all
				ps_rf = ps_rf_all
				
				ys = ys_all

				rfile = paste0("n", n, "_rho", rho*10, "_p0", p0, "_p", p, "_loop.RData")
				save(list=save_list, file=rfile)
			}else{
				rfile = paste0("n", n, "_rho", rho*10, "_p0", p0, "_p", p, "_loop.RData")
				load(rfile)
				ys = ys_test			
			}# end of if(merge)
			
			cur_rowname = paste0("n", n, "_rho", rho*10, "_p0", p0, "_p", p)
			y = t(apply(ys, 1, function(x){as.numeric(x)-1}))
			
			#########################
			## Proper Scoring Rules
			#########################

			index = 1:nrow(y)
				
			## CART
			ps = ps_cart
			scores_cart = t(sapply(index, function(i){colMeans(compute_scores(y[i,], ps[i,]))}))
			scores_mean_cart = colMeans(scores_cart)
			scores_sd_cart = apply(scores_cart, 2, sd)
			
			## RF
			ps = ps_rf
			scores_rf = t(sapply(index, function(i){colMeans(compute_scores(y[i,], ps[i,]))}))
			scores_mean_rf = colMeans(scores_rf)
			scores_sd_rf = apply(scores_rf, 2, sd)
			
			## CB
			ps = ps_CB
			scores_CB = t(sapply(index, function(i){colMeans(compute_scores(y[i,], ps[i,]))}))
			scores_mean_CB = colMeans(scores_CB)
			scores_sd_CB = apply(scores_CB, 2, sd)
			
			## logistic
			ps = ps_glm_only
			scores_logistic = t(sapply(index, function(i){colMeans(compute_scores(y[i,], ps[i,]))}))
			scores_mean_logistic = colMeans(scores_logistic)
			scores_sd_logistic = apply(scores_logistic, 2, sd)
			
			## l1-logistic
			ps = ps_la
			scores_l1_logistic = t(sapply(index, function(i){colMeans(compute_scores(y[i,], ps[i,]))}))
			scores_mean_l1_logistic = colMeans(scores_l1_logistic)
			scores_sd_l1_logistic = apply(scores_l1_logistic, 2, sd)
			
			## logistic refit
			ps = ps_glm
			scores_logistic_refit = t(sapply(index, function(i){colMeans(compute_scores(y[i,], ps[i,]))}))
			scores_mean_logistic_refit = colMeans(scores_logistic_refit)
			scores_sd_logistic_refit = apply(scores_logistic_refit, 2, sd)
			
			## fusion
			ps = ps_fusion
			scores_fusion = t(sapply(index, function(i){colMeans(compute_scores(y[i,], ps[i,]))}))
			scores_mean_fusion = colMeans(scores_fusion)
			scores_sd_fusion = apply(scores_fusion, 2, sd)
			
			
			scores_mean = rbind(scores_mean_CB, scores_mean_cart, scores_mean_rf,
						scores_mean_l1_logistic, scores_mean_logistic_refit, 
						scores_mean_fusion)
			scores_sd = rbind(scores_sd_CB, scores_sd_cart, scores_sd_rf,
						scores_sd_l1_logistic, scores_sd_logistic_refit, 
						scores_sd_fusion)
			
			# extract AUCs
			AUC = c(mean(index_CBs[,2]), mean(index_carts[,2]), mean(index_rfs[,2]),
				mean(index_las[,2]), mean(index_glms[,2]), mean(index_fusions[,2]))
			AUC_sd = c(sd(index_CBs[,2]), sd(index_carts[,2]), sd(index_rfs[,2]),
				sd(index_las[,2]), sd(index_glms[,2]), sd(index_fusions[,2]))
			
			rNames = c(paste0(cur_rowname, "_CB"), paste0(cur_rowname, "_CART"), paste0(cur_rowname, "_RF"),
					paste0(cur_rowname, "_l1_logistic"), paste0(cur_rowname, "_logistit_refit"),
					paste0(cur_rowname, "_FILTER"))
			scores_mean = cbind(AUC, scores_mean)
			scores_sd = cbind(AUC_sd, scores_sd)
			rownames(scores_mean) = rNames
			rownames(scores_sd) = rNames
			
			# format the output
			prec_out = 4
			out_matrix = matrix("", nrow=nrow(scores_mean), ncol=ncol(scores_mean))
			for(i in 1:nrow(scores_mean))
			{
				for(j in 1:ncol(scores_mean))
				{
					
					out_matrix[i,j] = paste0(formatC(round(scores_mean[i,j], prec_out), format='f', digits=prec_out))
					out_matrix[i,j] = paste0(formatC(round(scores_mean[i,j], prec_out), format='f', digits=prec_out), "(", formatC(round(scores_sd[i,j], prec_out), format='f', digits=prec_out), ")")
				}
			}
			rownames(out_matrix) = rNames
			colnames(out_matrix) = colnames(scores_mean)
			
			out_matrix_all = rbind(out_matrix_all, out_matrix)
		
		}## end of n

	}## end of rho
	
}## end of p0

# output the table
if(model=="thresh"){
	write.csv(out_matrix_all, file='Table2.csv')
}else if(model=="piece"){
	write.csv(out_matrix_all, file='Table3.csv')
}




