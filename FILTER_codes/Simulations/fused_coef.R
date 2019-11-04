
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

##read the parameters from the command line
args = commandArgs(T)
#args = c("100", "5", "500", "0", "1", "100", "100", "3", "0", "1")


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


##set the model parameters, in this code, K0 is useless
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




ps_fusion = NULL
index_fusions = NULL
rocs_fusion = list()
beta_fusions = NULL
beta0_fusions = NULL
Centers = list()
one_probs = NULL

save_list = c("ps_fusion", "index_fusions", "rocs_fusion", 
		"beta_fusions", "beta0_fusions", "Centers", "one_probs")

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


##start the main loop
start = proc.time()
for(ik in loop)
{

print(ik)

###############################################################################
##generate Xs
###############################################################################

##the Xs are generate from the iid standard norm distribution
set.seed(10000+ik)
X = matrix(rnorm(n*(p+q), mean=0, sd=0.5), nrow=n)


###############################################################################
##generate the cut-points
###############################################################################

class_number = 2
betas = c(0, rep(coef_level, class_number-1))

cut_quantiles = seq(1/class_number, by=1/class_number, length.out=class_number-1)
cut_points = qnorm(cut_quantiles)	
X_factor = data.frame(lapply(data.frame(X), function(x){cut(x, breaks=c(-Inf, cut_points, Inf), labels=1:class_number)}))
cnames = sapply(1:(p+q), function(x){paste("X", x, sep='')})

colnames(X_factor) = cnames
colnames(X) = cnames

##encode the X_factor to dummy variables
XX_true = dummy.data.frame(data=X_factor, names=cnames)
XX_int = XX_true


###############################################################################
##generate the coefficients 
###############################################################################

##add the first p0 true betas
beta_true = rep(betas, p0)

##add the 0s coefficients of the other p-p0 variables
beta_true = c(beta_true, rep(0, (p-p0)*class_number))

##add the second q0 true betas
beta_true = c(beta_true, rep(betas, q0))

##add the 0s coefficients of the other q-1 variables, each has 3 coefficients
beta_true = c(beta_true, rep(0, (q-q0)*class_number))


###############################################################################
##generate responses
###############################################################################

##generate the ys by the logistic model
linear_part = apply(XX_int, 1, function(x){sum(x*beta_true)})
intercept = -mean(linear_part)
p_y = intercept + linear_part
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

##set the parameters that CART needed
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
	dt = data.frame(X=X_b[,k], y_b)
	model = rpart(y_b~., method="class", data=dt, parms=list(split=crit), control=rpart.control(xval=cv, cp=complex))
	split = model$splits

	if(is.null(split))next
	split = split[which(split[,1]==nrow(X_b)),]

	var_split[[k]] = unique(c(var_split[[k]], as.numeric(split[4])))
}

##Bagging step with CART
for(b in 1:B)
{
	set.seed(b+1000*ik)
	index = sample(x=nrow(X_train), size=nrow(X_train), replace=TRUE)
	X_b = X_train[index,]
	y_b = y_train[index]
	rownames(X_b) = NULL

	for(i in 1:length(vars))
	{
		k = vars[i]
		dt = data.frame(X=X_b[,k], y_b)
		model = rpart(y_b~., method="class", data=dt, parms=list(split=crit), control=rpart.control(xval=cv, cp=complex))
		split = model$splits
	
		if(is.null(split))next
		split = split[which(split[,1]==nrow(X_b)),]
			
		var_split[[k]] = unique(c(var_split[[k]], as.numeric(split[4])))
	}
}


##take the mean of the cut point for each variable
centers = sapply(var_split, mean)
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


##record the index computed in this loop
index_fusions = rbind(index_fusions, index_fusion)
ps_fusion = rbind(ps_fusion, p_fusion)
roc_fusion = roc(p_fusion, y_test2)
rocs_fusion = c(rocs_fusion, list(roc_fusion))
one_probs = c(one_probs, one_prob)
beta_fusions = c(beta_fusions, list(beta_hat))
beta0_fusions = c(beta0_fusions, beta0_hat)


save(list=save_list, file=file_name)
}

end = proc.time()
end - start

warnings()
