

library(cluster)
library(dummies)
library(rpart)
library(glmnet)
library(AUC)
library(randomForest)



###############################################################################
##set the parameters
###############################################################################


main_path = "/home/FILTER_codes/Real_data/Threshold Effect CI/"
if(!dir.exists(main_path)){
	dir.create(main_path, recursive=TRUE)
}
setwd(main_path)




# M: number of experiments' replications
# Sep: divide the M replications into Sep parts, then we can do several independent experiments simultaneously
# round: the roundth part among Sep parts
# B: number of replications in the Bagging procedure in the FILTER
# K0: number of estimated threhold points for clustering after Bagging procedure in the FILTER
# method: method for constructing confidence interval for threshold point
# downsampling: logical, do the subsampling or not to construct the confidence interval
# cluster_type: clustering method after Bagging procedure in the FILTER
# CV: number of cross-validation folds for CART
# all_data: logical, use all data to construct the confidence interval or not

##the M loops will divided into Sep part, and this program will run the roundth part
#args = commandArgs(T)
#args = c("1", "1", "0", "100", "1", "subsampling", "TRUE")
args = c("900", "15", "1", "100", "1", "subsampling", "TRUE")

M = as.numeric(args[1])
Sep = as.numeric(args[2])
if(Sep==1){
	round = 0
}else{
	round = as.numeric(args[3])
}
B = as.numeric(args[4])
K0 = as.numeric(args[5])
method = args[6]
downsampling = as.logical(args[7])
cluster_type = "kmeans"
CV = 10
all_data = TRUE

##read the data
d = read.csv("real_data.csv")

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

sample_ratio = 1
if(downsampling){
	zeros = which(d$y==0)
	ones = which(d$y==1)

	set.seed(1)
	zeros_samples = sample(zeros, size=length(ones))
	set.seed(2)
	ones_samples = ones
	samples_index = sort(c(zeros_samples, ones_samples))

	data = d[samples_index,]
}else{
	data = d
}
n = nrow(data)

out_path = paste0(main_path, "CI_outputs/")
if(!dir.exists(out_path)){
	dir.create(out_path, recursive=TRUE)
}
setwd(out_path)

##the index of factor variables, which are no need to be discreted
factor_vars = c(2, 25)

if(round==0){
	loop = (1):(M)
}else{
	loop = (M/Sep * (round-1) + 1):(M/Sep * round)
}
if(method=="asymptotic"){
	loop = 1
}
print(loop)


savelist = c("Betas_l", "Betas_u", "Centers")

start = proc.time()
block_sizes = round(n^seq(from=0.95, to=0.75, length.out=10))
#block_sizes = n
for(m in block_sizes)
{

    file_name = paste0("Real_result_K", K0, "_B", B, "_", method, "_m", m, "_", round, ".RData")
    print(file_name)


    Centers = NULL
    Betas_l = NULL
    Betas_u = NULL

    ##repeat M in total resampling to evaluate the average the performance of the different methods

    for(ik in loop)
    {
        print(paste0("Block size: ", m, ", Loop: ", ik))
    
        if(method=="subsampling"){
            print("subsampling")
    
            if(all_data){
                train_index = 1:n
            }else{
                set.seed(1000+(ik-1)*2)
                train_index = sample(x=n, size=m, replace=F)
            }
    
            X = data.frame(data[,-1])
            y = as.factor(data$y)
            for(i in factor_vars)
            {
                X[,i] = as.factor(X[,i])
            }
    
            X_train_all = data.frame(X[train_index,])
            y_train_all = y[train_index]
        }else if(method=="asymptotic"){
    
            ## all data used for training, no test set
            X = data.frame(data[,-1])
            y = as.factor(data$y)
            for(i in factor_vars)
            {
                X[,i] = as.factor(X[,i])
            }
    
            X_train_all = data.frame(X)
            y_train_all = y
        }
    
    
        var_index = 1:ncol(X)
        one_prob = (table(y_train_all)/sum(table(y_train_all)))[2]
    
    
        ##data used for other methods
        X_train = X_train_all[,var_index]
        y_train = y_train_all
    
        ##############################################################
        ##Using CART+Bagging to get the cut points estimation
        ##############################################################
        if(is.null(factor_vars)){
            vars = c(1:ncol(X_train_all))
        }else{
            vars = c(1:ncol(X_train_all))[-factor_vars]
        }
    
        cv = CV
        complex = 0.0001
        crit = "gini"
    
        ##with some variables always in the model for estimation's robust
        all_vars = vars
        vars_remain = c(1,3,4,6)
        vars = setdiff(vars, vars_remain)
    
        X_train_all = X_train_all[,all_vars]
        y_train_all = y_train_all
    
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
    
        print("start CB")
        ##Bagging step with CART
        for(b in 1:B)
        {
            set.seed(b+100*ik)
            index = sample(x=nrow(X_train), size=nrow(X_train), replace=TRUE)
            X_b = X_train[index,]
            y_b = y_train[index]
            rownames(X_b) = NULL
    
            dt = data.frame(X_b, y_b)
            model = rpart(y_b~., method="class", data=dt, parms=list(split=crit), control=rpart.control(xval=cv, cp=complex))
    
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
    
    
        ##clustering the cut points by k-means
        #print("Start Clustering")
    
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
                set.seed(i*11+j+1111*ik)
                if(cluster_type=="kmeans"){
                        kc = kc + sort(as.vector(clara(data.frame(x), k=as.numeric(K[j]), metric="euclidean", samples=round(length(x)*0.8), rngR=TRUE)$medoids))
                }else if(cluster_type=="kmedians"){
                        kc = kc + sort(as.vector(clara(data.frame(x), k=as.numeric(K[j]), metric="manhattan", samples=round(length(x)*0.8), rngR=TRUE)$medoids))
                }
                }
                centers = c(centers, list(kc/I))
            }
        }
    
        Centers = rbind(Centers, unlist(centers))
    
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
    
        if(method=="asymptotic"){
            X_model = X_tilde
            y_model = y
        }else{
            X_model = X_tilde[train_index,]
            y_model = y[train_index]
        }
    
        ##############################################################
        ##Using penalized likelihood to estimate the coefficients
        ##with fusion penalty
        ##############################################################
    
        print("Start fusion")
        ##cross validation to choose best lambda
        nfold = 5
        tune = "auc"
    
        betas_l = NULL
        betas_u = NULL
        for(i in 1:ncol(X))
        {
            if(i %in% factor_vars)  next
    
            block = (2*i-1):(2*i)
            lambdas = exp(seq(log(0.001), log(5), length.out=100))
    
            X_model_single = X_model[,block]
            cv = cv.glmnet(x=X_model_single, y=y_model, family="binomial", type.measure="auc", nfolds=nfold, lambda=lambdas)
            lambda = cv$lambda.min
    
            ##using glmnet with lasso penalty to fit the model for transformed data
            model = glmnet(x=X_model_single, y=y_model, family="binomial", alpha=1, lambda=lambdas)
    
            gamma_hat = coef(model, s=lambda)
            beta_hat = as.vector(D[block, block]%*%gamma_hat[-1])
            names(beta_hat) = colnames(XXX)[block]
            beta0_hat = gamma_hat[1]
    
            beta_l = beta0_hat
            beta_u = beta0_hat + beta_hat[2]
    
            betas_l = c(betas_l, beta_l)
            betas_u = c(betas_u, beta_u)
        }
    
        ds = unlist(centers)
        names(betas_l) = names(ds)
        names(betas_u) = names(ds)
    
        Betas_l = rbind(Betas_l, betas_l)
        Betas_u = rbind(Betas_u, betas_u)
    
        save(list=savelist, file=file_name)
    }# end of ik
}# end of m

end = proc.time()
end - start


