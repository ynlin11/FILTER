

####################################################################################################
# most of the following functions are from the R package 'murphydiagram'
# we just need some intermediate results


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


# Input checker needed below
input.check <- function(f, y, t, alpha){
  if ( (length(t) > 1) | (length(alpha) > 1) | (length(f) != length(y)) | any(c(!is.vector(f), !is.vector(y), !is.vector(t), !is.vector(alpha))) ) stop("invalid input") 
}

# Extremal score for quantiles
S.quantile <- function(f, y, t, alpha){
  input.check(f, y, t, alpha)
  ((y < f) - alpha) * ((t < f) - (t < y))
}

# Take the positive part of a vector
pos <- function(x) x*(x >= 0)

# Extremal score for expectiles
S.expectile <- function(f, y, t, alpha){
  input.check(f, y, t, alpha)
  c1 <- abs((y < f) - alpha) 
  c2 <- pos(y-t) - pos(f-t) - (y-f)*(t<f)
  out <- c1 * c2
  out
}
  
# Newey-West VCV matrix
# u: vector of data
# prewhite: logical, should prewhitening be done?
# k: truncation lag for autocorrelations. If set to NULL, will be chosen automatically.
# meth: either "qs" (Quadratic Spectral, Andrews 1991) or anything else (Bartlett Kernel, Newey/West)
vHAC <- function(u, prewhite = FALSE, k = NULL, meth = "qs"){
  
  if (!is.matrix(u)) u <- as.matrix(u)
  
  n <- nrow(u)
  nreg <- ncol(u)
  rho <- sigma <- matrix(0, nreg, 1)
  
  # do a VAR(1) prewhitening
  if (prewhite == TRUE){
    reg.x <- matrix(u[1:(n-1), ], n-1, nreg)
    reg.y <- matrix(u[2:n, ], n-1, nreg)
    aux <- lm(reg.y~reg.x-1)
    beta <- matrix(unname(aux$coefficients), nreg, nreg)
    v <- matrix(unname(aux$residuals), n-1, nreg)
  } else {
    v <- u
    beta <- matrix(0, nreg, nreg)
  }
  nv <- nrow(v)
  
  # choose nr of lags (if not provided)
  if (is.null(k)){
    
    for (i in 1:nreg){
      aux <- lm(v[2:nv, i]~v[1:(nv-1), i]-1)
      rho[i] <- unname(aux$coefficients)
      sigma[i] <- sum(unname(aux$residuals)^2) / nv
    }
    
    if (meth == "qs"){
	  # See Eq. (6.4) on page 835 of Andrews (1991) -> Note that his sigma^2 corresponds to our sigma
      top <- sum( (4*(rho^2) * (sigma^2)) / ((1-rho)^8) )
      bot <- sum( (sigma^2) / ((1-rho)^4) )
      k <- ceiling(1.3221*((top/bot)*n)^(0.2))    
    } else {
      top <- sum( (4*(rho^2) * (sigma^2)) / (((1-rho)^6)*((1+rho)^2)) )
      bot <- sum( (sigma^2) / ((1-rho)^4) )
      k <- min(c(ceiling(1.1447*((top/bot)*n)^(1/3)), round(0.5*n)))
    }
    
  }
  
  # compute HAC
  vcv <- (t(v) %*% v) / (n-1)
  
  if (k > 0){
    if (meth == "qs"){
      del <- ((6 * pi)/(5 * k)) * (1:(n-1))
      w <- 3 * (sin(del) / del - cos(del)) / (del^2)  
      if (prewhite == FALSE){
        mlag <- n - 1 
      } else {
        mlag <- nv - 1
      }
    } else {
      w <- 1 - (1:k)/(k+1)
      mlag <- k
    }
    for (i in 1:mlag){
      cov <- t(v[(i+1):nv, , drop = FALSE]) %*% (v[1:(nv-i), , drop = FALSE]) / (n-1)
      vcv <- vcv + w[i]*(cov + t(cov))
    }
  }
  
  d <- solve(diag(nreg) - t(beta))
  hac <- d %*% vcv %*% t(d)
  
  return(list(hac = hac, k = k))  
  
}

get_grid <- function(f1, f2, y, functional, alpha, literal = TRUE, lhdiff = 1e-10){
  if (literal == TRUE){
    # Choose range of t's (see Corollary 2a, 2b in paper)
    if (functional == "quantile"){
      tseq <- sort(unique(c(f1, f2, y)))
    } else if (functional == "expectile") {
	  # NOTE: Use y's for Murphy diagrams also in case alpha = 1/2, even if not strictly needed
      aux1 <- sort(unique(c(f1, f2, y)))
      aux2 <- sort(unique(c(f1, f2))) - lhdiff
      tseq <- sort(unique(c(aux1, aux2)))   
    }
  } else {
    # Simply choose equally spaced grid of fixed length
    aux <- c(f1, f2, y)
#    tseq <- seq(from = min(aux) - 0.01*sd(aux), to = max(aux) + 0.01*sd(aux), length.out = 500)
    tseq <- seq(from = 0, to = 1, length.out = 500)
  }
  tseq
}


murphydiagram_diff2 <- function(f1, f2, y, functional = "expectile", alpha = 0.5, equally_spaced = FALSE, lag_truncate = 0, conf_level = 0.95){
  cex.gen <- 1.6
  # Some variables
  nobs <- length(y)
  scl <- abs(qnorm(0.5*(1-conf_level)))
  # Define function
  if (functional == "expectile"){
    g <- function(f, t) S.expectile(f, y, t, alpha)
  } else if (functional == "quantile"){
    g <- function(f, t) S.quantile(f, y, t, alpha)
  } else {
    stop("Please choose either expectile or quantile functional")
  }
  # Choose range of t's
  tseq <- get_grid(f1, f2, y, functional, alpha, literal = 1 - equally_spaced)
  df <- data.frame(tseq = tseq, s1 = numeric(length(tseq)), s2 = numeric(length(tseq)), lb = numeric(length(tseq)), ub = numeric(length(tseq)))
  for (j in 1:length(tseq)){
    aux1 <- g(f1, tseq[j])
    aux2 <- g(f2, tseq[j])
    df[j, 2:3] <- c(mean(aux1), mean(aux2))
    aux.v <- try(vHAC(aux1 - aux2, k = lag_truncate, meth = "bartlett")$hac/nobs)
    # HAC estimator won't invert for some t in the tails -> Set variance to zero in these cases
    if (class(aux.v) == "try-error") aux.v <- 0
    df[j, 4:5] <- mean(aux1-aux2) + c(-1, 1)*scl*sqrt(aux.v)
  }
    
  # Plot: Score difference + confidence interval
  if (all(y %in% c(0, 1))){
    xx <- c(-0.05, 1.05)
  } else {
    xx <- c(min(tseq) - 0.1, max(tseq) + 0.1)
  }

  return(df)
}


murphydiagram2 <- function(f1, f2, y, functional = "expectile", alpha = 0.5, labels = c("Method 1", "Method 2"), 
						  colors = NULL, equally_spaced = FALSE){
  cex.gen <- 1.6
  # Define function for extremal score
  if (functional == "expectile"){
    g <- function(f, t) S.expectile(f, y, t, alpha)
  } else if (functional == "quantile"){
    g <- function(f, t) S.quantile(f, y, t, alpha)
  } else {
    stop("Please choose either expectile or quantile functional")
  }
  # Unless specified otherwise: Use colors as in paper
  if (is.null(colors)) colors <- c("#D55E00", "#56B4E9", "#000000")
  # Grid of theta values
  tseq <- get_grid(f1, f2, y, functional, alpha, literal =  1 - equally_spaced)
  # Data frame with score entries for all theta values
  df <- data.frame(tseq = tseq, s1 = numeric(length(tseq)), s2 = numeric(length(tseq)))
  for (j in 1:length(tseq)){
    aux1 <- g(f1, tseq[j])
    aux2 <- g(f2, tseq[j])
    df[j, 2:3] <- c(mean(aux1), mean(aux2))
  }
  # Plot: Scores for both methods
  if (all(y %in% c(0, 1))){
    xx <- c(-0.05, 1.05)
  } else {
    xx <- c(min(tseq) - 0.1, max(tseq) + 0.1)
  }

  return(df)  
}


####################################################################################################




library(scoringRules)
library(murphydiagram)
#install.packages(c("scoringRules", "murphydiagram"))

#############################################################################
## Real Data
#############################################################################


main_path = "/home/FILTER_codes/Real_data/"
if(!dir.exists(main_path)){
	dir.create(main_path, recursive=TRUE)
}
setwd(main_path)

Tp = "insample" # 'outsample' can produce plots in Figure 3, 'insample' can produce Fig 4 and 5 in the Appendix
# the following parameters should match the model you fitted by real_data_prediction.R
tune_type = "optimal"
B = 100
round = 0
fold = 5
alpha = 0.9 # parameter in the Murphy diagram
cluster_type = "kmeans"
tune = "auc"
functional = "expectile"
K = 0 # actually, this parameter is useless, in the following, we will use the choice of K=6

name = paste0("K", K, "_B", B, "_", cluster_type, "_", tune, "_Folds", fold, "_", tune_type)
print(paste0(name, "_alpha", alpha*100))


setwd(main_path)
file = paste0("Real_result_", name, ".RData")
load(file)

# create the output directory
if(Tp=="insample"){
	out_path = paste0(main_path, "insamples/", functional, "/")
}else if(Tp=="outsample"){
	out_path = paste0(main_path, "outsamples/", functional, "/")
}

# if output directory already exists, clear it
if(!dir.exists(out_path)){
	dir.create(out_path, recursive=T)
}
setwd(out_path)

current_path = paste0(out_path, name, "/alpha", alpha*100, "/")
if(!dir.exists(current_path)){
	dir.create(current_path, recursive=T)
}
setwd(current_path)

if(Tp=="outsample"){
	y = t(apply(ys_out, 1, function(x){as.numeric(x)-1}))
}else if(Tp=="insample"){
	y = t(apply(ys_in, 1, function(x){as.numeric(x)-1}))
}


###########################################################################
## Murphy diagram, this will create plots in Fig 4 and 5 when Tp=INsample
###########################################################################
w = 12
h = 10

#### average experiments' results
index = 1:nrow(y)
cex.gen <- 1.6
colors <- c("#D55E00", "#56B4E9", "#000000")
xx <- c(-0.05, 1.05)


print("FILTER vs logistic")
## FILTER vs logistic
labels = c("Filter", "Logistic")
df = matrix(0, nrow=1000, ncol=3)
for(i in index)
{
	if(Tp=="outsample"){
		temp = murphydiagram2(f1=ps_fusion_out$'6'[i,], f2=ps_glm_only_out[i,], y=y[i,], functional=functional, alpha=alpha, labels=labels, equally_spaced=TRUE)
	}else if(Tp=="insample"){
		temp = murphydiagram2(f1=ps_fusion_in$'6'[i,], f2=ps_glm_only_in[i,], y=y[i,], functional=functional, alpha=alpha, labels=labels, equally_spaced=TRUE)
	}
	df = df + temp
}
df = df / length(index)

pdf("FILTER_vs_logistic(only)_murphy_ave.pdf", width=w, height=h)
matplot(x = df[,1], y = df[,2:3], type = "l", lty = 1, lwd = 4, xlab = "Parameter p", ylab = "", bty = "n", cex.lab = 2.2, 
          cex.axis = 2, xlim = xx, ylim = c(0, 1.2*max(df[,2:3])), col = colors)
abline(h = 0, lty = 2)
legend("top", labels, col = colors, lwd = 4, bty = "n", horiz = TRUE, cex = 2)
dev.off()

## CI
for(i in index)
{
	pdf(paste0("FILTER_vs_logistic(only)_murphy_diff_", i, ".pdf"), width=w, height=h)
	
	if(Tp=="outsample"){
		temp = murphydiagram_diff2(f1=ps_fusion_out$'6'[i,], f2=ps_glm_only_out[i,], y=y[i,], functional=functional, alpha=alpha, equally_spaced=TRUE)
	}else if(Tp=="insample"){
		temp = murphydiagram_diff2(f1=ps_fusion_in$'6'[i,], f2=ps_glm_only_in[i,], y=y[i,], functional=functional, alpha=alpha, equally_spaced=TRUE)
	}

	df = temp

	matplot(x = df$tseq, y = df[, 4:5], type = "n", ylab = "", xlab = "Parameter p",
	        bty = "n", cex.axis = 2, cex.lab = 2.2, xlim = xx, col = 1)
	polygon(c(df$t, rev(df$t)), c(df$ub, rev(df$lb)), col = "grey", border = NA)
	lines(x = df$t, y = (df$s1-df$s2), type = "l", col = 1, lwd = 2.5)
	abline(h = 0, lty = 2)
	dev.off()
}




print("FILTER vs l1-logistic")
## FILTER vs l1-logistic
labels = c("Filter", expression("l"[1]*"-logistic"))
df = matrix(0, nrow=1000, ncol=3)
for(i in index)
{
	if(Tp=="outsample"){
		temp = murphydiagram2(f1=ps_fusion_out$'6'[i,], f2=ps_la_out[i,], y=y[i,], functional=functional, alpha=alpha, labels=labels, equally_spaced=TRUE)
	}else if(Tp=="insample"){
		temp = murphydiagram2(f1=ps_fusion_in$'6'[i,], f2=ps_la_in[i,], y=y[i,], functional=functional, alpha=alpha, labels=labels, equally_spaced=TRUE)
	}
	df = df + temp
}
df = df / length(index)

pdf("FILTER_vs_logistic(l1)_murphy_ave.pdf", width=w, height=h)
matplot(x = df[,1], y = df[,2:3], type = "l", lty = 1, lwd = 4, xlab = "Parameter p", ylab = "", bty = "n", cex.lab = 2.2, 
          cex.axis = 2, xlim = xx, ylim = c(0, 1.2*max(df[,2:3])), col = colors)
abline(h = 0, lty = 2)
legend("top", labels, col = colors, lwd = 4, bty = "n", horiz = TRUE, cex = 2)
dev.off()

## CI

for(i in index)
{
	pdf(paste0("FILTER_vs_logistic(l1)_murphy_diff_", i, ".pdf"), width=w, height=h)

	if(Tp=="outsample"){
		temp = murphydiagram_diff2(f1=ps_fusion_out$'6'[i,], f2=ps_la_out[i,], y=y[i,], functional=functional, alpha=alpha, equally_spaced=TRUE)
	}else if(Tp=="insample"){
		temp = murphydiagram_diff2(f1=ps_fusion_in$'6'[i,], f2=ps_la_in[i,], y=y[i,], functional=functional, alpha=alpha, equally_spaced=TRUE)
	}
	df = temp

	matplot(x = df$tseq, y = df[, 4:5], type = "n", ylab = "", xlab = "Parameter p",
	        bty = "n", cex.axis = 2, cex.lab = 2.2, xlim = xx, col = 1)
	polygon(c(df$t, rev(df$t)), c(df$ub, rev(df$lb)), col = "grey", border = NA)
	lines(x = df$t, y = (df$s1-df$s2), type = "l", col = 1, lwd = 2.5)
	abline(h = 0, lty = 2)
	dev.off()
}



print("FILTER vs logistic-refit")
## FILTER vs logistic-refit
labels = c("Filter", "logistic-refit")
df = matrix(0, nrow=1000, ncol=3)
for(i in index)
{
	if(Tp=="outsample"){
		temp = murphydiagram2(f1=ps_fusion_out$'6'[i,], f2=ps_glm_out[i,], y=y[i,], functional=functional, alpha=alpha, labels=labels, equally_spaced=TRUE)
	}else if(Tp=="insample"){
		temp = murphydiagram2(f1=ps_fusion_in$'6'[i,], f2=ps_glm_in[i,], y=y[i,], functional=functional, alpha=alpha, labels=labels, equally_spaced=TRUE)
	}
	df = df + temp
}
df = df / length(index)

pdf("FILTER_vs_logistic_refit_murphy_ave.pdf", width=w, height=h)
matplot(x = df[,1], y = df[,2:3], type = "l", lty = 1, lwd = 4, xlab = "Parameter p", ylab = "", bty = "n", cex.lab = 2.2, 
          cex.axis = 2, xlim = xx, ylim = c(0, 1.2*max(df[,2:3])), col = colors)
abline(h = 0, lty = 2)
legend("top", labels, col = colors, lwd = 4, bty = "n", horiz = TRUE, cex = 2)
dev.off()

## CI
for(i in index)
{
	pdf(paste0("FILTER_vs_logistic_refit_murphy_diff_", i, ".pdf"), width=w, height=h)

	if(Tp=="outsample"){
		temp = murphydiagram_diff2(f1=ps_fusion_out$'6'[i,], f2=ps_glm_out[i,], y=y[i,], functional=functional, alpha=alpha, equally_spaced=TRUE)
	}else if(Tp=="insample"){
		temp = murphydiagram_diff2(f1=ps_fusion_in$'6'[i,], f2=ps_glm_in[i,], y=y[i,], functional=functional, alpha=alpha, equally_spaced=TRUE)
	}
	df = temp

	matplot(x = df$tseq, y = df[, 4:5], type = "n", ylab = "", xlab = "Parameter p",
	        bty = "n", cex.axis = 2, cex.lab = 2.2, xlim = xx, col = 1)
	polygon(c(df$t, rev(df$t)), c(df$ub, rev(df$lb)), col = "grey", border = NA)
	lines(x = df$t, y = (df$s1-df$s2), type = "l", col = 1, lwd = 2.5)
	abline(h = 0, lty = 2)
	dev.off()
}




print("FILTER vs CART")
## FILTER vs CART
labels = c("Filter", "CART")
df = matrix(0, nrow=1000, ncol=3)
for(i in index)
{
	if(Tp=="outsample"){
		temp = murphydiagram2(f1=ps_fusion_out$'6'[i,], f2=ps_cart_out[i,], y=y[i,], functional=functional, alpha=alpha, labels=labels, equally_spaced=TRUE)
	}else if(Tp=="insample"){
		temp = murphydiagram2(f1=ps_fusion_in$'6'[i,], f2=ps_cart_in[i,], y=y[i,], functional=functional, alpha=alpha, labels=labels, equally_spaced=TRUE)
	}
	df = df + temp
}
df = df / length(index)

pdf("FILTER_vs_cart_murphy_ave.pdf", width=w, height=h)
matplot(x = df[,1], y = df[,2:3], type = "l", lty = 1, lwd = 4, xlab = "Parameter p", ylab = "", bty = "n", cex.lab = 2.2, 
          cex.axis = 2, xlim = xx, ylim = c(0, 1.2*max(df[,2:3])), col = colors)
abline(h = 0, lty = 2)
legend("top", labels, col = colors, lwd = 4, bty = "n", horiz = TRUE, cex = 2)
dev.off()

## CI
for(i in index)
{
	pdf(paste0("FILTER_vs_cart_murphy_diff_", i, ".pdf"), width=w, height=h)
	
	if(Tp=="outsample"){
		temp = murphydiagram_diff2(f1=ps_fusion_out$'6'[i,], f2=ps_cart_out[i,], y=y[i,], functional=functional, alpha=alpha, equally_spaced=TRUE)
	}else if(Tp=="insample"){
		temp = murphydiagram_diff2(f1=ps_fusion_in$'6'[i,], f2=ps_cart_in[i,], y=y[i,], functional=functional, alpha=alpha, equally_spaced=TRUE)
	}
	df = temp

	matplot(x = df$tseq, y = df[, 4:5], type = "n", ylab = "", xlab = "Parameter p",
	        bty = "n", cex.axis = 2, cex.lab = 2.2, xlim = xx, col = 1)
	polygon(c(df$t, rev(df$t)), c(df$ub, rev(df$lb)), col = "grey", border = NA)
	lines(x = df$t, y = (df$s1-df$s2), type = "l", col = 1, lwd = 2.5)
	abline(h = 0, lty = 2)
	dev.off()
}


print("FILTER vs CB")
## FILTER vs CB
labels = c("Filter", "CB")
df = matrix(0, nrow=1000, ncol=3)
for(i in index)
{
	if(Tp=="outsample"){
		temp = murphydiagram2(f1=ps_fusion_out$'6'[i,], f2=ps_CB_out[i,], y=y[i,], functional=functional, alpha=alpha, labels=labels, equally_spaced=TRUE)
	}else if(Tp=="insample"){
		temp = murphydiagram2(f1=ps_fusion_in$'6'[i,], f2=ps_CB_in[i,], y=y[i,], functional=functional, alpha=alpha, labels=labels, equally_spaced=TRUE)
	}
	df = df + temp
}
df = df / length(index)

pdf("FILTER_vs_CB_murphy_ave.pdf", width=w, height=h)
matplot(x = df[,1], y = df[,2:3], type = "l", lty = 1, lwd = 4, xlab = "Parameter p", ylab = "", bty = "n", cex.lab = 2.2, 
          cex.axis = 2, xlim = xx, ylim = c(0, 1.2*max(df[,2:3])), col = colors)
abline(h = 0, lty = 2)
legend("top", labels, col = colors, lwd = 4, bty = "n", horiz = TRUE, cex = 2)
dev.off()

## CI

for(i in index)
{
	pdf(paste0("FILTER_vs_CB_murphy_diff_", i, ".pdf"), width=w, height=h)
	
	if(Tp=="outsample"){
		temp = murphydiagram_diff2(f1=ps_fusion_out$'6'[i,], f2=ps_CB_out[i,], y=y[i,], functional=functional, alpha=alpha, equally_spaced=TRUE)
	}else if(Tp=="insample"){
		temp = murphydiagram_diff2(f1=ps_fusion_in$'6'[i,], f2=ps_CB_in[i,], y=y[i,], functional=functional, alpha=alpha, equally_spaced=TRUE)
	}
	df = temp

	matplot(x = df$tseq, y = df[, 4:5], type = "n", ylab = "", xlab = "Parameter p",
	        bty = "n", cex.axis = 2, cex.lab = 2.2, xlim = xx, col = 1)
	polygon(c(df$t, rev(df$t)), c(df$ub, rev(df$lb)), col = "grey", border = NA)
	lines(x = df$t, y = (df$s1-df$s2), type = "l", col = 1, lwd = 2.5)
	abline(h = 0, lty = 2)
	dev.off()
}


print("FILTER vs RF")
## FILTER vs RF
labels = c("Filter", "RF")
df = matrix(0, nrow=1000, ncol=3)
for(i in index)
{
	if(Tp=="outsample"){
		temp = murphydiagram2(f1=ps_fusion_out$'6'[i,], f2=ps_rf_out[i,], y=y[i,], functional=functional, alpha=alpha, labels=labels, equally_spaced=TRUE)
	}else if(Tp=="insample"){
		temp = murphydiagram2(f1=ps_fusion_in$'6'[i,], f2=ps_rf_in[i,], y=y[i,], functional=functional, alpha=alpha, labels=labels, equally_spaced=TRUE)
	}
	df = df + temp
}
df = df / length(index)

pdf("FILTER_vs_rf_murphy_ave.pdf", width=w, height=h)
matplot(x = df[,1], y = df[,2:3], type = "l", lty = 1, lwd = 4, xlab = "Parameter p", ylab = "", bty = "n", cex.lab = 2.2, 
          cex.axis = 2, xlim = xx, ylim = c(0, 1.2*max(df[,2:3])), col = colors)
abline(h = 0, lty = 2)
legend("top", labels, col = colors, lwd = 4, bty = "n", horiz = TRUE, cex = 2)
dev.off()

## CI
for(i in index)
{
	pdf(paste0("FILTER_vs_rf_murphy_diff_", i, ".pdf"), width=w, height=h)

	if(Tp=="outsample"){
		temp = murphydiagram_diff2(f1=ps_fusion_out$'6'[i,], f2=ps_rf_out[i,], y=y[i,], functional=functional, alpha=alpha, equally_spaced=TRUE)
	}else if(Tp=="insample"){
		temp = murphydiagram_diff2(f1=ps_fusion_in$'6'[i,], f2=ps_rf_in[i,], y=y[i,], functional=functional, alpha=alpha, equally_spaced=TRUE)
	}
	df = temp

	matplot(x = df$tseq, y = df[, 4:5], type = "n", ylab = "", xlab = "Parameter p",
	        bty = "n", cex.axis = 2, cex.lab = 2.2, xlim = xx, col = 1)
	polygon(c(df$t, rev(df$t)), c(df$ub, rev(df$lb)), col = "grey", border = NA)
	lines(x = df$t, y = (df$s1-df$s2), type = "l", col = 1, lwd = 2.5)
	abline(h = 0, lty = 2)
	dev.off()
}


###############################################################################
## merge all curves, this will create plots in Fig 3 when Tp=outsample
###############################################################################

grid_length = 500

print("FILTER vs all")
## FILTER vs all
df = matrix(0, nrow=grid_length, ncol=3)
for(i in index)
{
	if(Tp=="outsample"){
		temp = murphydiagram2(f1=ps_fusion_out$'6'[i,], f2=ps_glm_out[i,], y=y[i,], functional=functional, alpha=alpha, equally_spaced=TRUE)
	}else if(Tp=="insample"){
		temp = murphydiagram2(f1=ps_fusion_in$'6'[i,], f2=ps_glm_in[i,], y=y[i,], functional=functional, alpha=alpha, equally_spaced=TRUE)
	}
	df = df + temp
}
df = df / length(index)

murphy_FILTER = df[,2]
murphy_logistic = df[,3]

df = matrix(0, nrow=grid_length, ncol=3)
for(i in index)
{
	if(Tp=="outsample"){
		temp = murphydiagram2(f1=ps_fusion_out$'6'[i,], f2=ps_glm_only_out[i,], y=y[i,], functional=functional, alpha=alpha, equally_spaced=TRUE)
	}else if(Tp=="insample"){
		temp = murphydiagram2(f1=ps_fusion_in$'6'[i,], f2=ps_glm_only_in[i,], y=y[i,], functional=functional, alpha=alpha, equally_spaced=TRUE)
	}
	df = df + temp
}
df = df / length(index)

murphy_logistic_refit = df[,3]

df = matrix(0, nrow=grid_length, ncol=3)
for(i in index)
{
	if(Tp=="outsample"){
		temp = murphydiagram2(f1=ps_fusion_out$'6'[i,], f2=ps_la_out[i,], y=y[i,], functional=functional, alpha=alpha, equally_spaced=TRUE)
	}else if(Tp=="insample"){
		temp = murphydiagram2(f1=ps_fusion_in$'6'[i,], f2=ps_la_in[i,], y=y[i,], functional=functional, alpha=alpha, equally_spaced=TRUE)
	}
	df = df + temp
}
df = df / length(index)

murphy_l1_logistic = df[,3]

df = matrix(0, nrow=grid_length, ncol=3)
for(i in index)
{
	if(Tp=="outsample"){
		temp = murphydiagram2(f1=ps_fusion_out$'6'[i,], f2=ps_cart_out[i,], y=y[i,], functional=functional, alpha=alpha, equally_spaced=TRUE)
	}else if(Tp=="insample"){
		temp = murphydiagram2(f1=ps_fusion_in$'6'[i,], f2=ps_cart_in[i,], y=y[i,], functional=functional, alpha=alpha, equally_spaced=TRUE)
	}
	df = df + temp
}
df = df / length(index)

murphy_cart = df[,3]

df = matrix(0, nrow=grid_length, ncol=3)
for(i in index)
{
	if(Tp=="outsample"){
		temp = murphydiagram2(f1=ps_fusion_out$'6'[i,], f2=ps_CB_out[i,], y=y[i,], functional=functional, alpha=alpha, equally_spaced=TRUE)
	}else if(Tp=="insample"){
		temp = murphydiagram2(f1=ps_fusion_in$'6'[i,], f2=ps_CB_in[i,], y=y[i,], functional=functional, alpha=alpha, equally_spaced=TRUE)
	}
	df = df + temp
}
df = df / length(index)

murphy_CB = df[,3]

df = matrix(0, nrow=grid_length, ncol=3)
for(i in index)
{
	if(Tp=="outsample"){
		temp = murphydiagram2(f1=ps_fusion_out$'6'[i,], f2=ps_rf_out[i,], y=y[i,], functional=functional, alpha=alpha, equally_spaced=TRUE)
	}else if(Tp=="insample"){
		temp = murphydiagram2(f1=ps_fusion_in$'6'[i,], f2=ps_rf_in[i,], y=y[i,], functional=functional, alpha=alpha, equally_spaced=TRUE)
	}
	df = df + temp
}
df = df / length(index)

murphy_rf = df[,3]
x = df[,1]

murphys = cbind(murphy_FILTER, murphy_l1_logistic, murphy_logistic_refit, murphy_cart, murphy_CB, murphy_rf)


pdf("murphy_merged.pdf", width=w, height=h)

par(mar=c(6,6,5,2))
cex.gen = 1.2
matplot(x=x, y=murphys, type = "l", lty=1:6, , lwd=2, xlab="", ylab="", bty="n", cex.lab=cex.gen, 
          cex.axis=cex.gen, xlim=xx, ylim=c(0, 1.2*max(murphys)), font.axis=1, cex.axis=2.2,
	    col=c("black", "red", "blue", "purple", "green", "magenta"))
mtext(side=1, text="Parameter p", line=3.5, cex=2.5)
mtext(side=2, text="Expected Score", line=3.5, cex=2.5)

abline(h = 0, lty = 2)
legend("topright", legend=c(paste0("Filter-6"), 
				expression("l"[1]*"-logistic"), "logistic-refit", "CART", "CB", "RF"),
	col=c("black", "red", "blue", "purple", "green", "magenta"), lty=1:6, cex=2, lwd=2)
dev.off()

best_matrix = matrix(0, ncol=6, nrow=grid_length)
best_index = apply(murphys, 1, which.min)
for(i in 1:grid_length)
{
	best_matrix[i, best_index[i]] = 1
}

pdf("best_forecast_ids.pdf", width=w, height=h)
par(mar=c(6,6,5,2))
image(x=x, y=1:6, z=best_matrix, xlab="", ylab="", 
	col=c("white", "dodgerblue3"), bty="n", main="", xaxt='n', yaxt='n')
axis(side=1, at=seq(0, 1, by=0.2), line=0.2, lwd=1, cex.axis=2.2)
axis(side=2, at=1:7, tick=FALSE, lwd=1, cex.axis=2.2)
mtext(side=1, text="Parameter p", line=3.5, cex=2.5)
mtext(side=2, text="Method", line=3.5, cex=2.5)
dev.off()



