
# set the working directory
main_path = "/home/FILTER_codes/Simulations/coefs/"
if(!dir.exists(main_path)){
	dir.create(main_path, recursive=TRUE)
}
setwd(main_path)

##############################################################################
## merge the results from fused_coef.R for different parts if needed
##############################################################################
merge = FALSE

if(merge){

    rho = 0
    p0 = 5
    p = 500
    files = dir(pattern="*.RData")

    ns = seq(200, 800, by=100)
    loop = 10

    for(i in 1:length(ns))
    {
        print(i)

        beta_fusions_all = NULL
        beta0_fusions_all = NULL
        Centers_all = NULL
        index_fusions_all = NULL
        one_probs_all = NULL
        ps_fusion_all = NULL
        rocs_fusion_all = NULL

        for(j in 1:loop)
        {
            file = files[loop*(i-1)+j]
            print(file)
            load(file)
            beta_fusions_all = c(beta_fusions_all, beta_fusions)
            beta0_fusions_all = c(beta0_fusions_all, beta0_fusions)
            Centers_all = c(Centers_all, Centers)
            index_fusions_all = rbind(index_fusions_all, index_fusions)
            one_probs_all = c(one_probs_all, one_probs)
            ps_fusion_all = rbind(ps_fusion_all, ps_fusion)
            rocs_fusion_all = c(rocs_fusion_all, rocs_fusion)

            file.remove(file)
        }

        beta_fusions = beta_fusions_all
        beta0_fusions = beta0_fusions_all
        Centers = Centers_all
        index_fusions = index_fusions_all
        one_probs = one_probs_all
        ps_fusion = ps_fusion_all
        rocs_fusion = rocs_fusion_all

        save_list = c("beta_fusions", "beta0_fusions", "Centers", "index_fusions", 
                "one_probs", "ps_fusion", "rocs_fusion")

        save_file = paste0("n", ns[i], "_rho", rho*10, "_p0", p0, "_p", p, "_loop.RData")
        save(list=save_list, file=save_file)
    }
}# end of if



##############################################################################
## output the indexes and plots in 4.1.1
##############################################################################

files = dir(pattern="*.RData")


# parameters should match the settings in the fused_coef.R
p = 500
ns = seq(200, 800, by=100)
row_names = ns
Centers_hat_sd = NULL
Centers_hat = NULL
betas_hat_true = NULL
betas_hat_non = NULL
beta0 = NULL

for(file in files)
{
	file_split = unlist(strsplit(file, '_'))
	n = as.numeric(gsub(pattern='n', replacement='', file_split[1]))
	rho = as.numeric(gsub(pattern='rho', replacement='', file_split[2]))
	p0 = as.numeric(gsub(pattern='p0', replacement='', file_split[3]))
	p = as.numeric(gsub(pattern='p', replacement='', file_split[4]))
	true_centers_index = c(1:p0, (p+1):(p+p0))
	true_betas_index = c(1:(2*p0), (2*p+1):(2*p+2*p0))
		
	load(file)

	# extract the cut points' estimations	
	centers_hat = sapply(Centers, function(x){x[true_centers_index]})
	Centers_hat = c(Centers_hat, list(apply(centers_hat, 1, function(x){quantile(abs(x), probs=c(0.025, 0.5, 0.975))}) ) )
	Centers_hat_sd = c(Centers_hat_sd, list(apply(centers_hat, 1, sd)))

	# extract the coefficients' estimations
	coefs = sapply(beta_fusions, function(x)x)
	
	coefs_evens = coefs[seq(from=2, to=nrow(coefs), by=2),]
	betas_hat_true = c(betas_hat_true, list(coefs_evens[true_centers_index, ]))
	betas_hat_non = c(betas_hat_non, list(coefs_evens[-true_centers_index, ]))
}

names(Centers_hat_sd) = row_names
names(Centers_hat) = row_names
names(betas_hat_true) = row_names
names(betas_hat_non) = row_names


# the output matrix for Table 1
out = NULL
############################################
## Performance of cut point esitimation
############################################

# average performance
center_means = colMeans(sapply(Centers_hat, function(x)x[2,]))
center_ups = colMeans(sapply(Centers_hat, function(x)x[1,]))
center_lows = colMeans(sapply(Centers_hat, function(x)x[3,]))

center_means_sd = colMeans(sapply(Centers_hat_sd, function(x)x))
centers = sapply(1:length(center_means), function(x){paste(round(center_means[x], 3), '(', round(center_means_sd[x], 3), ')', sep='')})
names(centers) = row_names

out = rbind(out, centers)

###########################
## average behavior
###########################

# plot the left part in Figure 2
lwd = 10

fit = lm(y~x-1, data=data.frame(y=log(center_means), x=log(ns)))
slope = coef(fit)

pdf(file="Mean.pdf", width=10, height=10)
pch_cex = 4
par(mar=c(6,6,4,4))
plot(x=log(ns), y=log(center_means), type='p', xlab="", ylab="", cex.axis=2,
	font.axis=2, ylim=c(-4, -1), pch=15, cex=pch_cex)
lines(x=log(ns), y=log(ns^(slope)), type='l', lwd=lwd)
mtext("Log(n)", side=1, line=3.5, cex=2.5, font=2)
mtext("Logarithm of absolute bias", side=2, line=3.5, cex=2.5, font=2)

lines(x=log(ns), y=log(ns^(-1/2)), lty=4, lwd=lwd)
lines(x=log(ns), y=log(ns^(-1/3)), lty=2, lwd=lwd)

dev.off()



# seperate variable performance
center_means = sapply(Centers_hat, function(x)x[2,])
center_ups = sapply(Centers_hat, function(x)x[1,])
center_lows = sapply(Centers_hat, function(x)x[3,])


# plot the individual variables' results
centers = rbind(center_ups, center_means, center_lows)

lwd = 10
type = 'p'
ns = seq(200, 800, by=100)
max_log = max(center_means)
min_log = min(center_means)
max_log = exp(-1.5)
min_log = exp(-3.5)


# plot the middle part in Figure 2
fit = lm(y~x-1, data=data.frame(y=log(as.vector(center_means[1:5,])), x=log(rep(ns, 5))))
slope = coef(fit)

pdf(file="F1.pdf", width=10, height=10)
par(mar=c(6,6,4,4))
pch_cex = 4
lwd_pch = 4
plot(x=log(ns), y=log(center_means[1,]), type='p', xlab="", ylab="", cex.axis=2,
	font.axis=2, ylim=c(-4, -1), pch=0, cex=pch_cex, lwd=lwd_pch)
lines(x=log(ns), y=log(ns^(slope)), type='l', lwd=lwd)
lines(x=log(ns), y=log(ns^(-1/2)), lty=4, lwd=lwd)
lines(x=log(ns), y=log(ns^(-1/3)), lty=2, lwd=lwd)

lines(x=log(ns), y=log(center_means[2,]), type=type, col='red', pch=1, cex=pch_cex, lwd=lwd_pch)
lines(x=log(ns), y=log(center_means[3,]), type=type, col='green3', pch=2, cex=pch_cex, lwd=lwd_pch)
lines(x=log(ns), y=log(center_means[4,]), type=type, col='blue', pch=3, cex=pch_cex, lwd=lwd_pch)
lines(x=log(ns), y=log(center_means[5,]), type=type, col='darkorange2', pch=4, cex=pch_cex, lwd=lwd_pch)
mtext("Log(n)", side=1, line=3.5, cex=2.5, font=2)
mtext("Logarithm of absolute bias", side=2, line=3.5, cex=2.5, font=2)


legend("top", legend=c(expression("X"["1"]),
			expression("X"["2"]), expression("X"["3"]), expression("X"["4"]), expression("X"["5"])), bty="o", 
	pch=c(0:4), col=c('black', 'red', 'green3', 'blue', 'darkorange2'), cex=3, horiz=TRUE,
	lty=c(rep(NA,5)), lwd=3, x.intersp=0.1, y.intersp=0.1, border="white", seg.len=1)

lines(x=log(rep(845, 2)), y=c(-6,0))
dev.off()




# plot the right part in Figure 2
fit = lm(y~x-1, data=data.frame(y=log(as.vector(center_means[6:10,])), x=log(rep(ns, 5))))
slope = coef(fit)

pdf(file="F2.pdf", width=10, height=10)
par(mar=c(6,6,4,4))
pch_cex = 4
lwd_pch = 4
plot(x=log(ns), y=log(center_means[6,]), type='p', xlab="", ylab="", cex.axis=2,
	font.axis=2, ylim=c(-4, -1), pch=0, cex=pch_cex, lwd=lwd_pch)
lines(x=log(ns), y=log(ns^(slope)), type='l', lwd=lwd)
lines(x=log(ns), y=log(ns^(-1/2)), lty=4, lwd=lwd)
lines(x=log(ns), y=log(ns^(-1/3)), lty=2, lwd=lwd)

lines(x=log(ns), y=log(center_means[7,]), type=type, col='red', pch=1, cex=pch_cex, lwd=lwd_pch)
lines(x=log(ns), y=log(center_means[8,]), type=type, col='green3', pch=2, cex=pch_cex, lwd=lwd_pch)
lines(x=log(ns), y=log(center_means[9,]), type=type, col='blue', pch=3, cex=pch_cex, lwd=lwd_pch)
lines(x=log(ns), y=log(center_means[10,]), type=type, col='darkorange2', pch=4, cex=pch_cex, lwd=lwd_pch)
mtext("Log(n)", side=1, line=3.5, cex=2.5, font=2)
mtext("Logarithm of absolute bias", side=2, line=3.5, cex=2.5, font=2)

legend("top", legend=c(expression("X"["6"]),
			expression("X"["7"]), expression("X"["8"]), expression("X"["9"]), expression("X"["10"])), bty="o", 
	pch=c(0:4), col=c('black', 'red', 'green3', 'blue', 'darkorange2'), cex=3, horiz=TRUE,
	lty=c(rep(NA,5)), lwd=3, x.intersp=0.1, y.intersp=0.1, border="white", seg.len=1)

lines(x=log(rep(845, 2)), y=c(-6,0))
dev.off()





############################################
## Performance of coefficients esitimation
############################################

# set true coefficient level
true_beta_level = 3
betas_true = sapply(betas_hat_true, function(x){ colSums((x-true_beta_level)^2) })
betas_non = sapply(betas_hat_non, function(x){ colSums(x^2) })


# compute RSEs and RMSEs
RSEs = sapply(1:nrow(betas_true), function(x){sqrt( (betas_true[x] + betas_non[x]) )} )
RMSEs = sapply(1:nrow(betas_true), function(x){sqrt( (betas_true[x] + betas_non[x]) / (2*p))} )

RSEs = matrix(0, nrow=nrow(betas_true), ncol=ncol(betas_true))
RMSEs = matrix(0, nrow=nrow(betas_true), ncol=ncol(betas_true))
for(i in 1:nrow(RSEs))
{
	for(j in 1:ncol(RSEs))
	{
		RSEs[i,j] = sqrt( betas_true[i,j] + betas_non[i,j] )
		RMSEs[i,j] = sqrt( (betas_true[i,j] + betas_non[i,j]) / (2*p) )
	}
}
colnames(RSEs) = colnames(RSEs)
colnames(RMSEs) = colnames(RMSEs)

RSE_mean = colMeans(RSEs)
RSE_sd = apply(RSEs, 2, sd)

RMSE_mean = colMeans(RMSEs)
RMSE_sd = apply(RMSEs, 2, sd)

out = rbind(out, RSE=paste(round(RSE_mean, 3), '(', round(RSE_sd, 3), ')', sep=''), 
		RMSE=paste(round(RMSE_mean, 3), '(', round(RMSE_sd, 3), ')', sep=''))


# compute sensitivity and specificity
TPs = sapply(betas_hat_true, function(x){apply(x, 2, function(t){length(which(t!=0))})})
FPs = sapply(betas_hat_non, function(x){apply(x, 2, function(t){length(which(t!=0))})})
TNs = sapply(betas_hat_non, function(x){apply(x, 2, function(t){length(which(t==0))})})
FNs = sapply(betas_hat_true, function(x){apply(x, 2, function(t){length(which(t==0))})})

sensitivity = colMeans( TPs / (TPs + FNs) )
specificity = colMeans( TNs / (TNs + FPs) )

sensitivity_sd = apply(TPs / (TPs + FNs), 2, sd)
specificity_sd = apply(TNs / (TNs + FPs), 2, sd)

sen = paste(round(sensitivity, 3), '(', round(sensitivity_sd, 3), ')', sep='') 
spe = paste(round(specificity, 3), '(', round(specificity_sd, 3), ')', sep='') 

out = rbind(out, sen, spe)

# output the Table 1
rownames(out) = c("MB_t", "RSE_vs", "RMSE_vs", "SEN_vs", "SPE_vs")
write.csv(out, "Table1.csv")

