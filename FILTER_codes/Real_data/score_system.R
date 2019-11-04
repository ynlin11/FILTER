
library(caret)
library(cluster)
library(dummies)
library(rpart)
library(glmnet)


###############################################################################
##set the parameters
###############################################################################

main_path = "/home/FILTER_codes/Real_data/"
if(!dir.exists(main_path)){
	dir.create(main_path, recursive=TRUE)
}
setwd(main_path)

#args = commandArgs(T)
args = c("100", "6", "kmeans", "auc")


# B: number of replications in the Bagging procedure in the FILTER
# K0: number of estimated threhold points for clustering after Bagging procedure in the FILTER
# cluster_type: clustering method after Bagging procedure in the FILTER
# tune: criteria for cross-validation when choosing lambda for FILTER

B = as.numeric(args[1])
K0 = as.numeric(args[2])
cluster_type = args[3]
tune = args[4]


##read the data
d = read.csv("real_data.csv", check.names=FALSE)

d_main = d[,c(3,4,7,8,38,39)]
vars_main = 1:6
y = d[,36]
d = d[,-c(1:8,36:39)]
y1 = y
y[y1>7.8] = 1
y[y1<=7.8] = 0
y = as.factor(y)
d = cbind(d_main, d, y)
var_names = colnames(d)
colnames(d)[1:(ncol(d)-1)] = paste0("x",1:(ncol(d)-1))
levels(d$x2) = c(0,1)
d$x25 = as.factor(d$x25)
levels(d$x2) = c(0,1)
d = d[complete.cases(d),]
d = d[,c(ncol(d),1:(ncol(d)-1))]

data = d
n = nrow(data)


##the index of factor variables, which are no need to be discreted
factor_vars = c(2, 25)

X = data.frame(data[,-1])
y = as.factor(data$y)
n = nrow(X)
for(i in factor_vars)
{
    X[,i] = as.factor(X[,i])
}

var_index = 1:ncol(X)
one_prob = (table(y)/sum(table(y)))[2]


##using all the data to get the score system
X_train = X
y_train = y


##############################################################
##Using CART+Bagging to get the cut points estimation
##############################################################
if(is.null(factor_vars)){
	vars = c(1:ncol(X_train))
}else{
	vars = c(1:ncol(X_train))[-factor_vars]
}

# parameters for CART
cv = 10 
complex = 0.0001
crit = "gini"


##with some variables always in the model for estimation's the stable solution
all_vars = vars
vars_remain = 1:6
vars = setdiff(vars, vars_remain)


y_b = y_train
X_b = X_train

names = colnames(X_b)
var_split = NULL
for(i in 1:length(all_vars))
{
	var_split = c(var_split, list(NULL))
}
names(var_split) = names[all_vars]


for(i in 1:length(vars))
{
	k = vars[i]
	dt = data.frame(X_b[,c(vars_remain,k)], y_b)
	model = rpart(y_b~., method="class", data=dt, parms=list(split=crit), control=rpart.control(xval=cv, cp=complex))
	split = model$splits

	# get the first layer result from CART's procedure
	if(is.null(split))next
	split = split[which(split[,1]==nrow(X_b)),]
	
	# record all the split points found
	for(j in 1:nrow(split))
	{
		index = which(names(var_split)==(rownames(split))[j])
		var_split[[index]] = unique(c(var_split[[index]], as.numeric(split[j,4])))
	}

}

##Bagging step with CART
for(b in 1:B)
{
	print(b)
	
	# generate Bagging sample
	set.seed(b)
	index = sample(x=nrow(X_train), size=nrow(X_train), replace=TRUE)
	X_b = X_train[index,]
	y_b = y_train[index]
	rownames(X_b) = NULL

	# due to the independency, find the split points seperately for each variable
	for(i in 1:length(vars))
	{
		k = vars[i]
		dt = data.frame(X_b[,c(vars_remain,k)], y_b)
		model = rpart(y_b~., method="class", data=dt, parms=list(split=crit), control=rpart.control(xval=cv, cp=complex))
		split = model$splits

		# get the first layer result from CART's procedure	
		if(is.null(split))next
		split = split[which(split[,1]==nrow(X_b)),]
		
		# record all the split points found
		for(j in 1:nrow(split))
		{
			index = which(names(var_split)==(rownames(split))[j])
			var_split[[index]] = unique(c(var_split[[index]], as.numeric(split[j,4])))
		}
	}
}



##clustering the split points by k-means
print("Start Clustering")

K = rep(K0, length(var_split))
centers = NULL
I = 100
removed_variables = NULL
for(j in 1:length(var_split))
{
    x = var_split[[j]]
    kc = rep(0, K[j])
    if(length(x)==0){
	  removed_variables = c(removed_variables, j)
	  next
    }else if(length(x)<=K[j]){
        centers = c(centers, list(x))
    }else{
        for(i in 1:I)
        {
		set.seed(i*11+j+1111)
		if(cluster_type=="kmeans"){
	            kc = kc + sort(as.vector(clara(data.frame(x), k=as.numeric(K[j]), metric="euclidean", samples=round(length(x)*0.8), rngR=TRUE)$medoids))
		}else if(cluster_type=="kmedians"){
	            kc = kc + sort(as.vector(clara(data.frame(x), k=as.numeric(K[j]), metric="manhattan", samples=round(length(x)*0.8), rngR=TRUE)$medoids))
		}
        }
        centers = c(centers, list(kc/I))
    }
}

# in case there are some variables are useless in the Bagging procedure
if(!is.null(removed_variables)){
	names(centers) = names(var_split)[-removed_variables]
}else{
	names(centers) = names(var_split)
}

##encode the X into dummy variables via estimated cut point for each variable
if(!is.null(removed_variables)){
	centers_index = all_vars[-removed_variables]
}else{
	centers_index = all_vars
}
XX = X[,var_index]
for(i in 1:length(centers_index))
{
	k = centers_index[i]
	x = XX[,k]
	class = rep(0, length(x))
	splits = centers[[i]]
	for(j in 1:length(splits))
	{
		class[which(x<=splits[j])] = class[which(x<=splits[j])] + 1
	}
	XX[,k] = as.factor(length(splits) + 1 - class)
}
if(!is.null(removed_variables)){
	XX = XX[,-all_vars[removed_variables]]
}
XXX = dummy.data.frame(data=data.frame(XX), names=colnames(XX))

# record number of levels for each used variable
variable_levels = NULL
for(i in 1:ncol(XX)){
    variable_levels = c(variable_levels, length(levels(XX[,i])))
}



################################################################################
# set the baseline and the levels of each variable, the default using the first
# level as the baseline
################################################################################

# generate transform matrix T
base = rep(1, ncol(XX))

D = matrix(0, sum(variable_levels), sum(variable_levels))
last_length = 0

for(i in 1:length(variable_levels))
{
	sub_matrix = matrix(0, nrow=variable_levels[i], ncol=variable_levels[i])

	for(j in 1:variable_levels[i])
	{
		row_number = j
		if(j==base[i]){
			sub_matrix[row_number, j] = 1
		}else if(j>base[i]){           
			sub_matrix[row_number, j-1] = -1
			sub_matrix[row_number, j] = 1
		}else if(j<base[i]){
			sub_matrix[row_number, j+1] = -1
			sub_matrix[row_number, j] = 1
		}
	}

	inv_submatrix = solve(sub_matrix)

	D[(last_length+1):(last_length+variable_levels[i]),(last_length+1):(last_length+variable_levels[i])] = inv_submatrix
	last_length = last_length + variable_levels[i]
}

X_tilde = as.matrix(XXX) %*% D
X_tilde = apply(X_tilde, 2, function(x){x-mean(x)})

X_model = X_tilde
y_model = y



##############################################################
##Using penalized likelihood to estimate the coefficients
##with fusion penalty
##############################################################

##cross validation to choose best lambda
fold = 5

##using AUC as the criteria
set.seed(159)
cv = cv.glmnet(x=X_model, y=y_model, family="binomial", type.measure="auc", nfolds=fold)

# choose lambda with given criteria
tune_type = "optimal"
least = which(cv$lambda==cv$lambda.min)
if(tune_type=="onese"){
	rule = cv$cvm[least] - cv$cvsd[least] / sqrt(fold)
	nonzeros = cv$nzero[which(cv$cvm>rule)]
	jump = which(names(cv$nzero)==names(nonzeros)[1])
	lambda = cv$lambda[jump]
}else if(tune_type=="twose"){
	rule = cv$cvm[least] - 2*cv$cvsd[least] / sqrt(fold)
	nonzeros = cv$nzero[which(cv$cvm>rule)]
	jump = which(names(cv$nzero)==names(nonzeros)[1])
	lambda = cv$lambda[jump]
}else if(tune_type=="threese"){
	rule = cv$cvm[least] - 3*cv$cvsd[least] / sqrt(fold)
	nonzeros = cv$nzero[which(cv$cvm>rule)]
	jump = which(names(cv$nzero)==names(nonzeros)[1])
	lambda = cv$lambda[jump]
}else if(tune_type=="threesigma"){
	rule = cv$cvm[least] - 3*cv$cvsd[least]
	nonzeros = cv$nzero[which(cv$cvm>rule)]
	jump = which(names(cv$nzero)==names(nonzeros)[1])
	lambda = cv$lambda[jump]
}else if(tune_type=="optimal"){
	lambda = cv$lambda.min
}

# fit the model
model = glmnet(x=X_model, y=y_model, family="binomial", alpha=1)
gamma_hat = coef(model, s=lambda)

# get the estimated coefficients
beta_hat = as.vector(D%*%gamma_hat[-1])
names(beta_hat) = colnames(XXX)
beta0_hat = gamma_hat[1]


##################################################
##transform to a score system
##################################################

##combine data with the same adjoint coefficient levels
combined_results = list()
start_tag = 1
with_centers = sort(c(factor_vars, as.numeric(sapply(names(centers), function(x){gsub(pattern='x', replacement='', x)}))))
for(i in 1:length(variable_levels))
{
	levels = variable_levels[i]
	end_tag = start_tag + levels - 1
	betas = beta_hat[start_tag:end_tag]

	# j is used for the case there are some variable without any split point by CART procedure
	j = with_centers[i]
	if(j %in% factor_vars){
		center = NULL
	}else{
		center = centers[[which(names(centers)==paste0("x", j))]]
	}

	combined_beta = betas[1]
	combined_center = NULL

	for(b in 2:length(betas))
	{
		if(betas[b]!=betas[b-1]){
			combined_beta = c(combined_beta, betas[b])
			if(!is.null(center)){
				combined_center = c(combined_center, center[b-1])
			}
		}
	}

	sign = ifelse(all(combined_beta[-1]<0), -1, 1)

	combined_beta = combined_beta + abs(min(combined_beta))

	result = list(beta=combined_beta, center=combined_center, sign=sign)
	combined_results = c(combined_results, list(result))

	start_tag = start_tag + levels
}
names(combined_results) = colnames(XX)

# find the selected variables
selected = list()
for(i in 1:length(combined_results))
{
	if(length(combined_results[[i]]$beta)!=1){
		selected = c(selected, list(combined_results[[i]]))
		names(selected)[length(selected)] = names(combined_results)[i]
	}
}
selected_index = sapply(names(selected), function(x){as.numeric(gsub("x", "", x))})
length(selected_index) 	# number of all selected variables

# find the variables with more than 3 levels
thresh = 3
selected_index = NULL
for(i in 1:length(combined_results))
{
	if(i >= 9){
		j = i+1
	}else{
		j = i
	}

	if(length(combined_results[[i]]$beta)>thresh) selected_index = c(selected_index, j)
}
length(selected_index)  # number of selected variables with more than 3 levels


# format the selected variables' results to compute the scores
for(i in 1:length(selected))
{
	combined_center = NULL
	center_cur = selected[[i]]$center
	for(b in 1:length(center_cur))
	{
		if(b==1){
			combined_center = c(combined_center, paste("<", round(center_cur[b],2), sep=''))
			if(length(center_cur)>1)
				combined_center = c(combined_center, paste(round(center_cur[b],2), "~", round(center_cur[b+1],2), sep=''))
		}else if(b<length(center_cur)){
			combined_center = c(combined_center, paste(round(center_cur[b],2), "~", round(center_cur[b+1],2), sep=''))
		}	
		if(b==length(center_cur)){
			combined_center = c(combined_center, paste(">=", round(center_cur[b],2), sep=''))
		}
	}
	selected[[i]] = c(selected[[i]], range=list(combined_center))
}


##transform to the scores
max_betas = NULL
max_betas = sapply(selected, function(x){max_betas = c(max_betas, max(x$beta))})
ratio = 100 / sum(max_betas)

# compute the scores
scores = NULL
for(i in 1:length(selected))
{
	score = round(selected[[i]]$beta * ratio, 2)
	selected[[i]] = c(selected[[i]], score=list(score))
	scores = c(scores, selected[[i]]$score)
}


##output the scores
file_name = paste("score_", length(selected), ".txt", sep="")
if(file.exists(file_name)){
	file.remove(file_name)
}

# format the scores output
for(i in 1:length(selected))
{
	item = selected[[i]]
	index = selected_index[i]

	colname = rep("", length(item$range))
	colname[1] = var_names[index]
	out = data.frame(colname, item$range, round(item$score,2))

	if(!file.exists(file_name)){
		write.table(t(c("Variable Name", "Range", "Score")), file=file_name, row.names=TRUE, col.names=FALSE)
		write.table(out, file=file_name, row.names=TRUE, col.names=FALSE, append=TRUE)
	}else{
		write.table(out, file=file_name, col.names=FALSE, append=TRUE)
	}
}


##summary the individuals' scores
Scores = X[,selected_index]
for(i in 1:ncol(Scores))
{
	x = Scores[,i]
	class = rep(0, length(x))
	splits = c(-Inf, selected[[i]]$center, Inf)
	cut_result = cut(x, breaks=splits, labels=FALSE)
	Scores[, i] = (selected[[i]]$score)[cut_result]
}
ind_scores = round(rowSums(Scores), 2)

summary(ind_scores)
cut_off = quantile(ind_scores, probs=1-c(0.1, 0.05, 0.01, 0.001, one_prob))
round(cut_off, 2)

















