
minimize_volatility = function(Is, span_m, sd_k)
{
	n = ncol(Is)
	
	## running mean of span m due to the stochastic approximation
	for(i in 1:n)
	{
		index = seq(from=i-span_m, to=i)
		index = index[which(index>0 & index<=n)]
	
		Is[,i] = rowMeans(as.matrix(Is[,index]))
	}
	
	## compute volatilities
	VI = NULL
	for(i in 1:n)
	{
		index = seq(from=i-sd_k, to=i)
		index = index[which(index>0 & index<=n)]		
	
		VI = c(VI, sum(apply(as.matrix(Is[,index]), 1, sd)))
	}
	
	selected_m = as.numeric(colnames(Is)[which.min(VI)])
	return(list(selected_m=selected_m, VI=VI))
}

main_path = "/home/FILTER_codes/Real_data/Threshold Effect CI/CI_output/"
if(!dir.exists(main_path)){
	dir.create(main_path, recursive=TRUE)
}
setwd(main_path)


# n: sample size
# Sep: divide the M replications into Sep parts, then we can do several independent experiments simultaneously
# round: the roundth part among Sep parts
# B: number of replications in the Bagging procedure in the FILTER
# K0: number of estimated threhold points for clustering after Bagging procedure in the FILTER
# method: method for constructing confidence interval for threshold point


# the following parameters should match the parameters in the 'threshold_effect_CI.R'
n = 328
block_sizes = round(n^seq(from=0.95, to=0.75, length.out=10))
Sep = 15
Rounds = 1:Sep
K0 = 1
B = 100
method = "subsampling"
rate = "cube"
alpha = 0.05

## estimation from the complete data
m = n
round = 0
file_name = paste0("Real_result_K", K0, "_B", B, "_", method, "_m", m, "_", round, ".RData")
load(file_name)
betas_l_hat = as.numeric(Betas_l)
betas_u_hat = as.numeric(Betas_u)
centers_hat = as.numeric(Centers)

p = length(betas_l_hat)
var_names = colnames(Betas_l)

CIs = list()
CIs_direct = list()
for(m in block_sizes)
{
	Betas_l_all = NULL
	Betas_u_all = NULL
	Centers_all = NULL

	for(round in Rounds)
	{
		file_name = paste0("Real_result_K", K0, "_B", B, "_", method, "_m", m, "_", round, ".RData")
		load(file_name)
	
		Betas_l_all = rbind(Betas_l_all, Betas_l)
		Betas_u_all = rbind(Betas_u_all, Betas_u)
		Centers_all = rbind(Centers_all, Centers)
	}

	## use quantiles directly
	CI_l_direct = sapply(1:p, function(i){quantile(Betas_l_all[,i], probs=c(alpha/2, 0.5, 1-alpha/2))})
	colnames(CI_l_direct) = var_names
	CI_u_direct = sapply(1:p, function(i){quantile(Betas_u_all[,i], probs=c(alpha/2, 0.5, 1-alpha/2))})
	colnames(CI_u_direct) = var_names
	CI_d_direct = sapply(1:p, function(i){quantile(Centers_all[,i], probs=c(alpha/2, 0.5, 1-alpha/2))})
	colnames(CI_d_direct) = var_names
	current_CIs_direct = list(CI_l=CI_l_direct, CI_u=CI_u_direct, CI_d=CI_d_direct)
	CIs_direct = c(CIs_direct, list(current_CIs_direct))

	## set rate
	if(rate=="root"){
		tau_n = n^(1/2)
		tau_m = m^(1/2)
	}else if(rate=="cube"){
		tau_n = n^(1/3)
		tau_m = m^(1/3)
	}

	## solve CI for three parameters
	quantiles_l = sapply(1:p, function(i){quantile(tau_m*(Betas_l_all[,i]-betas_l_hat[i]), probs=c(alpha/2, 0.5, 1-alpha/2))})
	colnames(quantiles_l) = var_names
	CI_l = betas_l_hat - tau_n^(-1) * quantiles_l

	quantiles_u = sapply(1:p, function(i){quantile(tau_m*(Betas_u_all[,i]-betas_u_hat[i]), probs=c(alpha/2, 0.5, 1-alpha/2))})
	colnames(quantiles_u) = var_names
	CI_u = betas_u_hat - tau_n^(-1) * quantiles_u
	
	quantiles_d = sapply(1:p, function(i){quantile(tau_m*(Centers_all[,i]-centers_hat[i]), probs=c(alpha/2, 0.5, 1-alpha/2))})
	colnames(quantiles_d) = var_names
	CI_d = centers_hat - tau_n^(-1) * quantiles_d
	
	current_CIs = list(CI_l=CI_l, CI_u=CI_u, CI_d=CI_d)
	CIs = c(CIs, list(current_CIs))
}
names(CIs) = block_sizes
names(CIs_direct) = block_sizes

span_m = 3
sd_k = 3
ms_left = NULL
CIs_left = NULL
ms_right = NULL
CIs_right = NULL
ms_center = NULL
CIs_center = NULL

CIs_used = CIs_direct
for(i in 1:p)
{
	left = sapply(CIs_used, function(x){x$CI_l[,i]})
	m_star = minimize_volatility(left, span_m, sd_k)$selected_m
	ms_left = c(ms_left, m_star)
	CIs_left = cbind(CIs_left, left[,which(colnames(left)==m_star)])

	right = sapply(CIs_used, function(x){x$CI_u[,i]})
	m_star = minimize_volatility(right, span_m, sd_k)$selected_m
	ms_right = c(ms_right, m_star)
	CIs_right = cbind(CIs_right, right[,which(colnames(right)==m_star)])

	center = sapply(CIs_used, function(x){x$CI_d[,i]})
	m_star = minimize_volatility(center, span_m, sd_k)$selected_m
	ms_center = c(ms_center, m_star)
	CIs_center = cbind(CIs_center, center[,which(colnames(center)==m_star)])
}
colnames(CIs_left) = var_names
colnames(CIs_right) = var_names
colnames(CIs_center) = var_names


library(plotrix)

factor_index = c(2, 25)
var_names_all = read.csv("var_names.csv", colClasses="character")
var_names = var_names_all$Abr[-factor_index]

threshed_index = NULL
for(i in 1:p)
{
	if(max(CIs_left[c(1,3),i]) < min(CIs_right[c(1,3),i]) | 
		min(CIs_left[c(1,3),i]) > max(CIs_right[c(1,3),i])){
		threshed_index = c(threshed_index, i)
	}
}

sorted_index = c(threshed_index, (1:p)[-threshed_index])

ulim = max(rbind(CIs_left, CIs_right))
llim = min(rbind(CIs_left, CIs_right))


# plot the Fig 1 in the paper
w = 18
h = 10
pdf("CI comparison.pdf", width=w, height=h)

range = 1.05
par(mar=c(10,3,2,2))
CIs_content = CIs_left
cols = c(rep("blue", length(threshed_index)), rep("grey", p-length(threshed_index)))
lwds = 2
plotCI(x=1:p, y=CIs_content[2,sorted_index], ui=CIs_content[3,sorted_index], li=CIs_content[1,sorted_index],
	col=cols, ylim=c(llim*range, ulim*range), lwd=lwds, xlab="", cex.lab=2, cex=1.5, slty=2, gap=0, sfrac=0.006,
	ylab="", axes=F)
title(xlab="Variable", line=8, cex.lab=2)

par(new=T)
CIs_content = CIs_right
cols = c(rep("red", length(threshed_index)), rep("grey", p-length(threshed_index)))
lwds = 3
plotCI(x=1:p, y=CIs_content[2,sorted_index], ui=CIs_content[3,sorted_index], li=CIs_content[1,sorted_index],
	col=cols, lwd=lwds, xlab="", ylab="", ylim=c(llim*range, ulim*range), axes=F, cex=1.5, slty=1, gap=0, sfrac=0.006,
	pch=4)

x_labels = var_names[sorted_index]
axis(side=1, at=1:p, labels=rep("", p), tck=.015, lwd.ticks=1, lwd=1)
#axis(1, at=1:p, labels=x_labels, line=0, tick=F, cex.axis=1.5, las=2, adj=1)
text(1:p, par("usr")[3]-0.5, labels=x_labels, srt=45, pos=1, xpd=TRUE, adj=1, cex=1.5)
axis(side=2, at=c(-3:4), labels=c(-3:4), tck=.015, lwd.ticks=1, cex.axis=1.5, lwd=1)

abline(h=-2.605, lwd=1)
lines(x=c(0, 13.5), y=c(0,0), col="blue", lwd=2)
lines(x=c(13.5, 43), y=c(0,0), col="grey", lwd=2)
#abline(h=0, col="black", lwd=1)

cols = c("red", "blue")
#cols = "black"
legend("top", legend=c("Right Region", "Left Region"), pch=c(4, 1), 
	lwd=c(3, 2), col=cols, lty=c(1, 2), cex=1.5, x.intersp=0.1, seg.len=2,
	bty='n', horiz=TRUE)
dev.off()



## plot variables with threshold effect

threshed_index = NULL
for(i in 1:p)
{
	if(max(CIs_left[c(1,3),i]) < min(CIs_right[c(1,3),i]) | 
		min(CIs_left[c(1,3),i]) > max(CIs_right[c(1,3),i])){
		threshed_index = c(threshed_index, i)
	}
}


ulim = max(rbind(CIs_left[,threshed_index], CIs_right[,threshed_index]))
llim = min(rbind(CIs_left[,threshed_index], CIs_right[,threshed_index]))

w = 14
h = 10
pdf("CI comparison_threhold effect.pdf", width=w, height=h)

range = 1.05
par(mar=c(10,3,2,2))
CIs_content = CIs_left[,threshed_index]
plotCI(x=1:ncol(CIs_content), y=CIs_content[2,], ui=CIs_content[3,], li=CIs_content[1,],
	col="blue", ylim=c(llim*range, ulim*range), lwd=2, xlab="", cex.lab=2, cex=1.5, 
	ylab="", axes=F)
title(xlab="Variable", line=8, cex.lab=2)

par(new=T)
CIs_content = CIs_right[,threshed_index]
plotCI(x=1:ncol(CIs_content), y=CIs_content[2,], ui=CIs_content[3,], li=CIs_content[1,],
	col="red", lwd=3, xlab="", ylab="", ylim=c(llim*range, ulim*range), axes=F, cex=1.5,
	pch=4)

#x_labels = (1:43)[-factor_index]
x_labels = var_names[threshed_index]
axis(side=1, at=1:ncol(CIs_content), labels=rep("", ncol(CIs_content)), tck=.015, lwd.ticks=1, lwd=1)
#axis(1, at=1:ncol(CIs_content), labels=x_labels, line=0, tick=F, cex.axis=1.5, las=2)
text(1:ncol(CIs_content), par("usr")[3]-0.5, labels=x_labels, srt=45, pos=1, xpd=TRUE, adj=1, cex=1.5)
axis(side=2, at=c(-3:4), labels=c(-3:4), tck=.015, lwd.ticks=1, cex.axis=1.5, lwd=1)

abline(h=-2.605, lwd=1)
abline(h=0, col="blue", lwd=2)

legend("top", legend=c("Right Region", "Left Region"), pch=c(4, 1), 
	lwd=c(3, 2), col=c("red", "blue"), lty=0, cex=1.5, x.intersp=0.1, 
	bty='n', horiz=TRUE)
dev.off()




