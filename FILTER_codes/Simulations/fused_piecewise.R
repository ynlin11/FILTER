
library(MASS)
library(dummies)
library(rpart)
library(glmnet)
library(AUC)
library(randomForest)


##evaluate function for prediction performance
eva = function(y_test, y_pred, p_pred)
{
	library(AUC)

	table = table(y_test, y_pred)
	if(ncol(table)==1){
		if(colnames(table)=="0"){
			table = cbind(table, c(0, 0))
		}else if(colnames(table)=="1"){
			table = cbind(c(0, 0), table)
		}
	}
	acc = sum(diag(table)) / length(y_test)

	##Recall, sensitivity
	Sensitivity = table[2,2] / (table[2,2] + table[2,1])
		
	##Precision
	if((table[2,2] + table[1,2])==0){
		Precision = 0
	}else{
		Precision = table[2,2] / (table[2,2] + table[1,2])
	}
		
	##FDR, equals to (1 - precision)
	FDR = 1 - Precision
	
	##True negative rate, specificity
	Specificity = table[1,1] / (table[1,1] + table[1,2])
		
	##False Positive Rate, equals to (1 - specificity)
	FPR = table[1,2] / (table[1,2] + table[1,1])
	
	##AUC by ROC
	AUC = auc(roc(p_pred, labels=y_test))
	
	index = c(ACC=acc, AUC=AUC, Sensitivity=Sensitivity, Precision=Precision, FDR=FDR, Specificity=Specificity, FPR=FPR)
	
	return(index)
}




###############################################################################
##set the parameters
###############################################################################

##set the main path to restore the results
main_path = "/home/FILTER_codes/Simulations/"
if(!dir.exists(main_path)){
	dir.create(main_path, recursive=TRUE)
}
setwd(main_path)


args = commandArgs(T)
#args = c("50", "2", "10", "0", "6", "10", "3", "5", "0", "1")


# M: number of experiments' replications
# Sep: divide the M replications into Sep parts, then we can do several independent experiments simultaneously
# round: the roundth part among Sep parts
# B: number of replications in the Bagging procedure in the FILTER
# K0: number of estimated threhold points for clustering after Bagging procedure in the FILTER
# CV: number of cross-validation folds for CART
# n_train: half of the sample size
# n: sample size
# p0, q0: number of non-zeros in the true model
# p, q: number of all coefficients in the model
# rho: parameter for the covariance matrix of the noise vector
# coef_level: magnitude for true coefficients


##set the model parameters
n_train = as.numeric(args[1])
n = n_train*2
p0 = as.numeric(args[2])
q0 = p0
p = as.numeric(args[3])
q = p
rho = as.numeric(args[4])
K0 = as.numeric(args[5])
M = as.numeric(args[6])
B = as.numeric(args[7])
coef_level = as.numeric(args[8])
CV = 10

##the M loops will divided into Sep part, and this program will run the roundth part
round = as.numeric(args[9])
Sep = as.numeric(args[10])




ps_cart = NULL
index_carts = NULL
threshes_cart = NULL
rocs_cart = list()

ps_rf = NULL
index_rfs = NULL
threshes_rf = NULL
rocs_rf = list()

ps_CB = NULL
index_CBs = NULL
threshes_CB = NULL
rocs_CB = list()

ps_fusion = NULL
index_fusions = NULL
rocs_fusion = list()

ps_la = NULL
index_las = NULL
rocs_la = list()

ps_glm = NULL
index_glms = NULL
rocs_glm = list()

ps_glm_only = NULL
index_glms_only = NULL
rocs_glm_only = list()

Centers = list()
one_probs = NULL

save_list = c("ps_cart", "index_carts", "rocs_cart", 
		"ps_rf", "index_rfs", "rocs_rf",
		"ps_glm_only", "index_glms_only", "rocs_glm_only",
		"ps_CB", "index_CBs", "rocs_CB", "one_probs")
temp = c("ps_fusion", "index_fusions", "rocs_fusion", 
		"ps_la", "index_las", "rocs_la",
		"ps_glm", "index_glms", "rocs_glm")
save_list = c(save_list, temp, "Centers", "ys_test")


if(round==0){
	loop = (1):(M)
	file_name = paste("n", n, "_rho", rho*10, "_p0", p0, "_p", p, "_loop.RData", sep='')
}else{
	loop = (M/Sep * (round-1) + 1):(M/Sep * round)
	file_name = paste("n", n, "_rho", rho*10, "_p0", p0, "_p", p, "_loop", round, ".RData", sep='')
}


print(loop)
print(file_name)

##main effects, generated from a multivariates normal distribution,
##which has mu0=0 mean and (Sigma0)_{i,j}=\rho^{|i-j|} covariance matrix
mu0 = rep(0, p)

Sigma0 = diag(rep(0.5,p))
for(i in 1:(p-1))
{
	for(j in (i+1):p)
	{
		Sigma0[i,j] = rho^(abs(i-j))
	}
}
Sigma0 = Sigma0 + t(Sigma0)


start = proc.time()
for(ik in loop)
{

print(ik)

###############################################################################
##generate Xs
###############################################################################

set.seed(100+ik)
X = mvrnorm(n=n, mu0, Sigma0)

for(i in 1:q)
{
	set.seed(10000+ik*100+i)
	X = cbind(X, rd=rnorm(n))
}



###############################################################################
##generate the coefficients 
###############################################################################

class_number = 4

cut_quantiles = seq(1/class_number, by=1/class_number, length.out=class_number-1)
cut_points = qnorm(cut_quantiles)	

beta_true_forth = rep(0, p+q)
beta_true_left = rep(0, p+q)
beta_true_right = rep(0, p+q)
beta_true_middle = rep(0, p+q)

true_index = c(1:p0, (p+1):(p+q0))
beta_true_forth[true_index] = 5
beta_true_middle[true_index] = 10
beta_true_right[true_index] = 5


###############################################################################
##generate responses
###############################################################################

##generate the ys by the logistic model
if(class_number==3){
	model_sum = apply(X, 1, function(x){
		sum(sin(pi*x)*ifelse(x<cut_points[1] | x==cut_points[1], 1, 0)*beta_true_left + 
		(x-0.5)^2*ifelse(x>cut_points[2] | x==cut_points[2], 1, 0)*beta_true_right +
		x*ifelse(cut_points[1]<x & x<cut_points[2], 1, 0)*beta_true_middle )})
	intercept = -mean(model_sum)
	p_y = intercept + model_sum
}else if(class_number==2){
	model_sum = apply(X, 1, function(x){
		sum(sin(pi*x)*ifelse(x<cut_points[1] | x==cut_points[1], 1, 0)*beta_true_left + 
		(x-0.5)^2*ifelse(x>cut_points[1], 1, 0)*beta_true_right )})
	intercept = -mean(model_sum)
	p_y = intercept + model_sum
}else if(class_number==4){
	model_sum = apply(X, 1, function(x){
		sum(x*ifelse(x<cut_points[1] | x==cut_points[1], 1, 0)*beta_true_left + 
		sin(pi*x)*ifelse(cut_points[1]<x & x<cut_points[2], 1, 0)*beta_true_forth +
		(x-0.5)^2*ifelse(cut_points[2]<x & x<cut_points[3], 1, 0)*beta_true_middle + 
		x*ifelse(x>cut_points[3] | x==cut_points[3], 1, 0)*beta_true_right)})
	intercept = -mean(model_sum)
	p_y = intercept + model_sum
}


p_y = exp(p_y) / (1+exp(p_y))

set.seed(100000+ik*111)
y = as.factor(sapply(p_y, function(x){rbinom(1, size=1, prob=x)}))


####################################
##Start training
####################################

set.seed(1000+(ik-1)*7)
train_index = sample(x=n, size=n*0.8, replace=F)

X_train = data.frame(X[train_index,])
y_train = y[train_index]
X_test = data.frame(X[-train_index,])
y_test = y[-train_index]
ys_test = y_test

##recode the positive samples proportion
one_prob = (table(y_train)/sum(table(y_train)))[2]



##############################################################
##Using CART+Bagging to get the cut points estimation
##############################################################

##the index of factor variables, which are no need to be discreted
factor_vars = NULL
if(is.null(factor_vars)){
	vars = c(1:ncol(X_train))
}else{
	vars = c(1:ncol(X_train))[-factor_vars]
}

all_vars = vars
vars_remain = 1:3
vars = setdiff(vars, vars_remain)

cv = CV
complex = 0.00001
crit = "gini"


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

	if(is.null(split))next
	split = split[which(split[,1]==nrow(X_b)),]
	
	for(j in 1:nrow(split))
	{
		index = which(names(var_split)==(rownames(split))[j])
		var_split[[index]] = unique(c(var_split[[index]], as.numeric(split[j,4])))
	}

}


pred_prob = matrix(0, nrow=nrow(X_test), ncol=B)
y_t = as.numeric(y_test) - 1

##Bagging step with CART
for(b in 1:B)
{
	set.seed(b+1000*ik)
	index = sample(x=nrow(X_train), size=nrow(X_train), replace=TRUE)
	X_b = X_train[index,]
	y_b = y_train[index]
	rownames(X_b) = NULL

	model = rpart(y_b~., method="class", data=data.frame(X_b), parms=list(split=crit), control=rpart.control(xval=cv, cp=complex))
	pred_prob[,b] = predict(model, newdata=data.frame(X_test), type="prob")[,2]

	for(i in 1:length(vars))
	{
		k = vars[i]
		dt = data.frame(X_b[,c(vars_remain,k)], y_b)
		model = rpart(y_b~., method="class", data=dt, parms=list(split=crit), control=rpart.control(xval=cv, cp=complex))
		split = model$splits
	
		if(is.null(split))next
		split = split[which(split[,1]==nrow(X_b)),]
		
		for(j in 1:nrow(split))
		{
			index = which(names(var_split)==(rownames(split))[j])
			var_split[[index]] = unique(c(var_split[[index]], as.numeric(split[j,4])))
		}
	}
}


##record the result of CART+Bagging for prediction in this loop
p_CB = rowMeans(pred_prob)
roc_CB = roc(p_CB, y_test)
thresh_CB = roc_CB$cutoff[which.max(abs(roc_CB$fpr - roc_CB$tpr))]
y_CB = as.factor(ifelse(p_CB>thresh_CB, 1, 0))
index_CB = eva(y_test, y_CB, p_CB)

ps_CB = rbind(ps_CB, as.vector(p_CB))
index_CBs = rbind(index_CBs, index_CB)
rocs_CB = c(rocs_CB, list(roc_CB))
threshes_CB = c(threshes_CB, thresh_CB)



##clustering the cut points by k-means
K = rep(K0, length(var_split))

centers = NULL
I = 200
for(j in 1:length(var_split))
{
    x = var_split[[j]]
    kc = rep(0, K[j])
    if(length(x)<=K[j]){
        centers = c(centers, list(x))
    }else{
        for(i in 1:I)
        {
		set.seed(i*111+j+ik*1111)
            kc = kc + sort(kmeans(x, as.numeric(K[j]), algorithm="MacQueen", iter.max=25)$centers)
        }
        centers = c(centers, list(kc/I))
    }
}
names(centers) = names(var_split)
Centers = c(Centers, list(centers))

##encode the X into dummy variables via estimated cut point for each variable
XX = data.frame(X[,all_vars])

for(i in 1:length(all_vars))
{
    k = all_vars[i]
    x = XX[,k]
    class = rep(0, length(x))
    splits = centers[[i]]
    for(j in 1:length(splits))
    {
        class[which(x<=splits[j])] = class[which(x<=splits[j])] + 1
    }
    XX[,k] = as.factor(length(splits) + 1 - class)
}

XXX = dummy.data.frame(data=data.frame(XX), names=colnames(XX))

variable_levels = NULL
for(i in 1:ncol(XX)){
    variable_levels = c(variable_levels, length(levels(XX[,i])))
}


################################################################################
##set the baseline and the levels of each variable, the default using the first
##level as the baseline
################################################################################
base = rep(1, ncol(XX))
base = c(base[1], cumsum((variable_levels))[-length(variable_levels)]+base[-1])

D = matrix(0, sum(variable_levels), sum(variable_levels))
last_length = 0

for(i in 1:length(variable_levels))
{
    sub_matrix = matrix(0, nrow=variable_levels[i], ncol=variable_levels[i])
    
    if(i==1){
        base_index = base[i]
    }else{
        base_index = base[i] - sum(variable_levels[1:(i-1)]+1)
    }   

    if(base_index==1){
        for(j in 1:variable_levels[i])
        {
    
            row_number = j
            if(j==1){
                sub_matrix[row_number, j] = 1
                next
            }
            
            sub_matrix[row_number, j-1] = -1
            sub_matrix[row_number, j] = 1
        }
    }else if(base_index==(variable_levels[i]+1)){
        for(j in 1:variable_levels[i])
        {
    
            row_number = j
            if(j==variable_levels[i]){
                sub_matrix[row_number, j] = 1
                next
            }
            
            sub_matrix[row_number, j] = -1
            sub_matrix[row_number, j+1] = 1
        }       
    }else{
        for(j in 1:variable_levels[i])
        {
    
            row_number = j
            if(j==base_index | j==(base_index-1)){
                sub_matrix[row_number, j] = 1
                next
            }
            
            sub_matrix[row_number, j-1] = -1
            sub_matrix[row_number, j] = 1
        }       
    }
    
    inv_submatrix = solve(sub_matrix)

    D[(last_length+1):(last_length+variable_levels[i]),(last_length+1):(last_length+variable_levels[i])] = inv_submatrix
    last_length = last_length + variable_levels[i]
}

X_tilde = as.matrix(XXX) %*% D
X_tilde = apply(X_tilde, 2, function(x){x-mean(x)})


X_model = X_tilde[train_index,]
y_model = y[train_index]

X_test2 = X_tilde[-train_index, ]
y_test2 = y[-train_index]


##############################################################
##Using penalized likelihood to estimate the coefficients
##with fusion penalty
##############################################################

##cross validation to choose best lambda
fold = 5

##using AUC as the criteria
cv = cv.glmnet(x=X_model, y=y_model, family="binomial", type.measure="auc", nfolds=fold)

##using likelihood as the criteria
#cv = cv.glmnet(x=X_model, y=y_model, family="binomial", type.measure="deviance", nfolds=fold)

##using accuracy as the criteria
#cv = cv.glmnet(x=X_model, y=y_model, family="binomial", type.measure="class", nfolds=fold)



lambda = cv$lambda.min

##using glmnet with lasso penalty to fit the model for transformed data
model = glmnet(x=X_model, y=y_model, family="binomial", alpha=1)
gamma_hat = coef(model, s=lambda)

##beta_hat is the estimated coefficients
beta_hat = as.vector(D%*%gamma_hat[-1])
names(beta_hat) = colnames(XXX)

beta0_hat = gamma_hat[1]

##compute the index for prediction performance
p_fusion = as.vector(predict(model, newx=X_test2, type="response", s=lambda))
y_fusion = as.factor(ifelse(p_fusion>one_prob, 1, 0))
index_fusion = eva(y_test2, y_fusion, p_fusion)


##record the index for CART+Fusion method computed in this loop
index_fusions = rbind(index_fusions, index_fusion)
ps_fusion = rbind(ps_fusion, p_fusion)
roc_fusion = roc(p_fusion, y_test2)
rocs_fusion = c(rocs_fusion, list(roc_fusion))

one_probs = c(one_probs, one_prob)


##############################################################
##Other methods' performance
##############################################################

cv = 10 
complex = 0.00001
crit = "gini"


X_train = X[train_index,]
y_train = y[train_index]
X_test = X[-train_index,]
y_test = y[-train_index]


##CART

dt = data.frame(X_train, y_train)
model_cart = rpart(y_train~., method="class", data=dt, parms=list(split=crit), control=rpart.control(xval=cv, cp=complex))

p_cart = predict(model_cart, newdata=data.frame(X_test), type="prob")[,2]
roc_cart = roc(p_cart, y_test)
thresh_cart = roc_cart$cutoff[which.max(abs(roc_cart$fpr - roc_cart$tpr))]
y_cart = as.factor(ifelse(p_cart>thresh_cart, 1, 0))
index_cart = eva(y_test, y_cart, p_cart)

ps_cart = rbind(ps_cart, as.vector(p_cart))
index_carts = rbind(index_carts, index_cart)
rocs_cart = c(rocs_cart, list(roc_cart))
threshes_cart = c(threshes_cart, thresh_cart)


##RF
dt = data.frame(X_train, y_train)
prior = c(1-one_prob, one_prob)
names(prior) = c("0", "1")
model_rf = randomForest(y_train~., data=dt, ntree=500, cutoff=prior, xtest=data.frame(X_test), ytest=y_test)

y_rf = model_rf$test$predicted
p_rf = model_rf$test$votes[,2]
index_rf = eva(y_test, y_rf, p_rf)
roc_rf = roc(p_rf, y_test)

ps_rf = rbind(ps_rf, as.vector(p_rf))
index_rfs = rbind(index_rfs, index_rf)
rocs_rf = c(rocs_rf, list(roc_rf))


###Logistic(only)

glm_data = data.frame(X_train, y=y_train)
glm_model = glm(y~., data=glm_data, family="binomial")

p_glm_only= predict(glm_model, newdata=data.frame(X_test), type="response")
roc_glm_only = roc(p_glm_only, y_test)
y_glm_only = as.factor(ifelse(p_glm_only>one_prob, 1, 0))
index_glm_only = eva(y_test, y_glm_only, p_glm_only)

ps_glm_only = rbind(ps_glm_only, as.vector(p_glm_only))
index_glms_only = rbind(index_glms_only, index_glm_only)
rocs_glm_only = c(rocs_glm_only, list(roc_glm_only))



##Logistic(l1)
cv = cv.glmnet(x=X_train, y=y_train, family="binomial", type.measure="auc", nfolds=10)
model = glmnet(x=X_train, y=y_train, family="binomial", alpha=1, lambda=cv$lambda.min)

beta = as.numeric(model$beta)

p_la = predict(model, newx=X_test, type="response")
y_la = as.factor(ifelse(p_la>one_prob, 1, 0))
index_la = eva(y_test, y_la, p_la)
roc_la = roc(p_la, y_test)

ps_la = rbind(ps_la, as.vector(p_la))
index_las = rbind(index_las, index_la)
rocs_la = c(rocs_la, list(roc_la))


##Logistic(refit)
glm_data = data.frame(X=X_train[,which(beta!=0)], y=y_train)
glm_model = glm(y~., data=glm_data, family="binomial")

p_glm= predict(glm_model, newdata=data.frame(X=X_test[,which(beta!=0)]), type="response")
y_glm = as.factor(ifelse(p_glm>one_prob, 1, 0))
index_glm = eva(y_test, y_glm, p_glm)
roc_glm = roc(p_glm, y_test)

ps_glm = rbind(ps_glm, as.vector(p_glm))
index_glms = rbind(index_glms, index_glm)
rocs_glm = c(rocs_glm, list(roc_glm))


save(list=save_list, file=file_name)
}

end = proc.time()
end - start


##############################################################
##Compute the average performance for different methods
##############################################################

mean_CB=  colMeans(index_CBs)
sd_CB = apply(index_CBs, 2, sd)

mean_cart =  colMeans(index_carts)
sd_cart = apply(index_carts, 2, sd)

mean_rf =  colMeans(index_rfs)
sd_rf = apply(index_rfs, 2, sd)

mean_fusion =  colMeans(index_fusions)
sd_fusion = apply(index_fusions, 2, sd)

mean_la =  colMeans(index_las)
sd_la = apply(index_las, 2, sd)

mean_glm =  colMeans(index_glms)
sd_glm = apply(index_glms, 2, sd)

mean_glm_only =  colMeans(index_glms_only)
sd_glm_only = apply(index_glms_only, 2, sd)


result = rbind(mean_CB, sd_CB, mean_cart, sd_cart, mean_rf, sd_rf, 
		mean_fusion, sd_fusion, mean_la, sd_la, 
		mean_glm, sd_glm, mean_glm_only, sd_glm_only)

CB = NULL
CART = NULL
RF = NULL
fusion = NULL
lasso = NULL
glm = NULL
glm_only = NULL

for(i in 1:length(mean_CB))
{
	CB = c(CB, paste(round(mean_CB[i], 3), "(", round(sd_CB[i], 3), ")", sep=''))
	CART = c(CART, paste(round(mean_cart[i], 3), "(", round(sd_cart[i], 3), ")", sep=''))
	RF = c(RF, paste(round(mean_rf[i], 3), "(", round(sd_rf[i], 3), ")", sep=''))
	fusion = c(fusion, paste(round(mean_fusion[i], 3), "(", round(sd_fusion[i], 3), ")", sep=''))
	lasso = c(lasso, paste(round(mean_la[i], 3), "(", round(sd_la[i], 3), ")", sep=''))
	glm = c(glm, paste(round(mean_glm[i], 3), "(", round(sd_glm[i], 3), ")", sep=''))
	glm_only = c(glm_only, paste(round(mean_glm_only[i], 3), "(", round(sd_glm_only[i], 3), ")", sep=''))
}


result = rbind(CB, CART, RF, fusion, lasso, glm, glm_only)
colnames(result) = names(mean_CB)


#Results
file_name = paste("n", n, "_rho", rho*10, "_p0", p0, "_p", p, "_round", round, ".csv", sep='')
write.csv(result, file_name, row.names=TRUE)


file_name = paste("n", n, "_rho", rho*10, "_p0", p0, "_p", p, "_round", round, ".RData", sep='')
save(list=save_list, file=file_name)

