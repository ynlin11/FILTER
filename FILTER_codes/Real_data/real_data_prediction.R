
library(caret)
library(cluster)
library(dummies)
library(rpart)
library(glmnet)
library(AUC)
library(randomForest)


# evaluate function for prediction performance
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
	Precision = table[2,2] / (table[2,2] + table[1,2])
		
	##FDR, equals to (1 - precision)
	FDR = table[1,2] / (table[1,2] + table[2,2])	
	
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


main_path = "/home/FILTER_codes/Real_data/"
if(!dir.exists(main_path)){
	dir.create(main_path, recursive=TRUE)
}
setwd(main_path)

out_path = paste0(main_path, "outputs/")


##the M loops will divided into Sep part, and this program will run the roundth part
args = commandArgs(T)
#args = c("100", "0", "1", "100", "0", "kmeans", "auc", "TRUE", "5", "optimal")


# M: number of experiments' replications
# Sep: divide the M replications into Sep parts, then we can do several independent experiments simultaneously
# round: the roundth part among Sep parts
# B: number of replications in the Bagging procedure in the FILTER
# K0s: numbers of estimated threhold points for clustering after Bagging procedure in the FILTER
# CV: number of cross-validation folds for CART
# cluster_type: clustering method after Bagging procedure in the FILTER
# tune: criteria for cross-validation when choosing lambda for FILTER
# cv_train: logical, indicate whether do the cross-validation analysis for the data or not
# fold: number of croos-validation folds for the whole analysis procedure
# tune_type: tuning criteria to choose lambda for FILTER


M = as.numeric(args[1])
round = as.numeric(args[2])
Sep = as.numeric(args[3])
B = as.numeric(args[4])
K0 = as.numeric(args[5]) # actually, K0 is useless, since we will permute all the values in K0s as K0, this is just for debug
cluster_type = args[6]
tune = args[7]
cv_train = as.logical(args[8])
fold = as.numeric(args[9])
tune_type = args[10]
if(!(tune_type %in% c("optimal", "onese", "twose", "threese", "threesigma"))) stop("Error parameter in tune_type!")
file_name = paste0("Real_result_K", K0, "_B", B, "_", cluster_type, "_", tune, "_Folds", fold, "_", tune_type, ".RData")
K0s = c(2, 4, 6, 8)

CV = 10
print(file_name)

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

data = d
n = nrow(data)


out_path = paste0(main_path, "outputs/")
if(!dir.exists(out_path)){
	dir.create(out_path, recursive=TRUE)
}
setwd(out_path)



##the index of factor variables, which are no need to be discreted
factor_vars = c(2, 25)

one_probs = NULL

index_fusions = list(NULL, NULL, NULL, NULL)
ps_fusion_out = list(NULL, NULL, NULL, NULL)
ps_fusion_in = list(NULL, NULL, NULL, NULL)
rocs_fusion = list(list(), list(), list(), list())
threshes_fusion = list(NULL, NULL, NULL, NULL)
names(index_fusions) = K0s
names(ps_fusion_out) = K0s
names(ps_fusion_in) = K0s
names(rocs_fusion) = K0s
names(threshes_fusion) = K0s


ps_la_out = NULL
ps_la_in = NULL
index_las = NULL
threshes_la = NULL
rocs_la = list()

ps_glm_out = NULL
ps_glm_in = NULL
index_glms = NULL
threshes_glm = NULL
rocs_glm = list()

ps_glm_only_out = NULL
ps_glm_only_in = NULL
index_glms_only = NULL
threshes_glm_only = NULL
rocs_glm_only = list()

ps_cart_out = NULL
ps_cart_in = NULL
index_carts = NULL
threshes_cart = NULL
rocs_cart = list()

ps_rf_out = NULL
ps_rf_in = NULL
index_rfs = NULL
threshes_rf = NULL
rocs_rf = list()

ps_CB_out = NULL
ps_CB_in = NULL
index_CBs = NULL
threshes_CB = NULL
rocs_CB = list()

ys_out = NULL
ys_in = NULL


if(round==0){
	loop = (1):(M)
}else{
	loop = (M/Sep * (round-1) + 1):(M/Sep * round)
}

print(loop)


savelist = c("ik", "ps_fusion_out", "ps_fusion_in", "index_fusions", "rocs_fusion",
		"ps_la_in", "ps_la_out", "index_las", "rocs_la", "ps_glm_out", "ps_glm_in", "index_glms", "rocs_glm",
		"ps_glm_only_out", "ps_glm_only_in", "index_glms_only", "rocs_glm_only", "ys_in", "ys_out",
		"ps_cart_out", "ps_cart_in", "index_carts", "rocs_cart", "ps_rf_out", "ps_rf_in", "index_rfs", "rocs_rf",
		"ps_CB_out", "ps_CB_in", "index_CBs", "rocs_CB")


## create fold id
set.seed(fold)
flds = createFolds(d$y, k=fold, list=TRUE, returnTrain=FALSE)


##repeat M in total resampling to evaluate the average the performance of the different methods
start = proc.time()
for(ik in loop)
{
	print(paste0("Loop: ", ik))

	if(cv_train){
		# use cross-validation to evaluate methods
		print("CV_train")	

		cv_loop = (1:fold)[-ik]
		train_index = sort(unlist(sapply(cv_loop, function(i){flds[[i]]})))

		X = data.frame(data[,-1])
		y = as.factor(data$y)
		for(i in factor_vars)
		{
			X[,i] = as.factor(X[,i])
		}

		set.seed(111+ik)
		X_train_all = data.frame(X[train_index,])
		y_train_all = y[train_index]
		
		X_test = data.frame(X[-train_index,])
		y_test = y[-train_index]
	}else{
		##divide the data into training and testing data
		set.seed(1000+(ik-1)*2)
		train_index = sample(x=n, size=n*0.8, replace=F)
		
		X = data.frame(data[,-1])
		y = as.factor(data$y)
		for(i in factor_vars)
		{
			X[,i] = as.factor(X[,i])
		}
		
		X_train_all = data.frame(X[train_index,])
		y_train_all = y[train_index]
		
		X_test = data.frame(X[-train_index,])
		y_test = y[-train_index]
	}


	var_index = 1:ncol(X)
	one_prob = (table(y_train_all)/sum(table(y_train_all)))[2]


	##data used for other methods
	X_train = X_train_all[,var_index]
	y_train = y_train_all
	X_test = X_test[, var_index]


	##############################################################
	##Using CART+Bagging to get the cut points estimation
	##############################################################
	if(is.null(factor_vars)){
		vars = c(1:ncol(X_train))
	}else{
		vars = c(1:ncol(X_train))[-factor_vars]
	}
	
	# parameters for CART
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


	pred_probs = NULL
	fit_probs = NULL
	y_t = as.numeric(y_test) - 1

	print("start CB")
	##Bagging step with CART
	for(b in 1:B)
	{

		# generate Bagging sample
		set.seed(b+100*ik)
		index = sample(x=nrow(X_train), size=nrow(X_train), replace=TRUE)
		X_b = X_train[index,]
		y_b = y_train[index]
		rownames(X_b) = NULL

		dt = data.frame(X_b, y_b)
		model = rpart(y_b~., method="class", data=dt, parms=list(split=crit), control=rpart.control(xval=cv, cp=complex))

		fit_probs = cbind(fit_probs, predict(model, newdata=X_b, type="prob")[,2])
		pred_probs = cbind(pred_probs, predict(model, newdata=data.frame(X_test), type="prob")[,2])

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


	##record the result of CART+Bagging for prediction in this loop
	p_CB = rowMeans(pred_probs)
	p_CB_in = rowMeans(fit_probs)
	roc_CB = roc(p_CB, y_test)
	thresh_CB = one_prob
	y_CB = as.factor(ifelse(p_CB>thresh_CB, 1, 0))
	index_CB = eva(y_test, y_CB, p_CB)


	ps_CB_out = rbind(ps_CB_out, as.vector(p_CB))
	ps_CB_in = rbind(ps_CB_in, as.vector(p_CB_in))
	index_CBs = rbind(index_CBs, index_CB)
	rocs_CB = c(rocs_CB, list(roc_CB))
	threshes_CB = c(threshes_CB, thresh_CB)


	##clustering the cut points by k-means
	print("Start Clustering")


	K0_it = 1
	for(K0 in K0s)
	{
		print(paste0("K0: ", K0))
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

		variable_levels = NULL
		for(i in 1:ncol(XX)){
			variable_levels = c(variable_levels, length(levels(XX[,i])))
		}



		################################################################################
		##set the baseline and the levels of each variable, the default using the first
		##level as the baseline
		################################################################################
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

		X_model = X_tilde[train_index,]
		y_model = y[train_index]

		X_test2 = X_tilde[-train_index, ]
		y_test2 = y[-train_index]

		##############################################################
		##Using penalized likelihood to estimate the coefficients
		##with fusion penalty
		##############################################################

		print("Start fusion CV")
		##cross validation to choose best lambda
		nfold = 5

		if(tune=="auc"){
			##using AUC as the criteria
			cv = cv.glmnet(x=X_model, y=y_model, family="binomial", type.measure="auc", nfolds=nfold)
		}else if(tune=="likelihood"){
			##using likelihood as the criteria
			cv = cv.glmnet(x=X_model, y=y_model, family="binomial", type.measure="deviance", nfolds=nfold)
		}else if(tune=="acc"){
			##using accuracy as the criteria
			cv = cv.glmnet(x=X_model, y=y_model, family="binomial", type.measure="class", nfolds=nfold)
		}

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


		##using glmnet with lasso penalty to fit the model for transformed data
		model = glmnet(x=X_model, y=y_model, family="binomial", alpha=1)

		p_fusion = as.vector(predict(model, newx=X_test2, type="response", s=lambda))
		fit_fusion = as.vector(predict(model, newx=X_model, type="response", s=lambda))
		roc_fusion = roc(p_fusion, y_test2)
		thresh_fusion = one_prob
		y_fusion = as.factor(ifelse(p_fusion>thresh_fusion, 1, 0))
		index_fusion = eva(y_test2, y_fusion, p_fusion)

		index_fusions[[K0_it]] = rbind(index_fusions[[K0_it]], index_fusion)
		ps_fusion_out[[K0_it]] = rbind(ps_fusion_out[[K0_it]], p_fusion)
		ps_fusion_in[[K0_it]] = rbind(ps_fusion_in[[K0_it]], fit_fusion)
		rocs_fusion[[K0_it]] = c(rocs_fusion[[K0_it]], list(roc_fusion))
		threshes_fusion[[K0_it]] = c(threshes_fusion[[K0_it]], thresh_fusion)

		K0_it = K0_it + 1
	}## end of K0


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

	print("Start Others")
	##CART
	dt = data.frame(X_train, y_train)
	model_cart = rpart(y_train~., method="class", data=dt, parms=list(split=crit), control=rpart.control(xval=cv, cp=complex))

	p_cart = predict(model_cart, newdata=data.frame(X_test), type="prob")[,2]
	fit_cart = predict(model_cart, newdata=data.frame(X_train), type="prob")[,2]
	roc_cart = roc(p_cart, y_test)
	thresh_cart = one_prob
	y_cart = as.factor(ifelse(p_cart>thresh_cart, 1, 0))
	index_cart = eva(y_test, y_cart, p_cart)

	ps_cart_out = rbind(ps_cart_out, as.vector(p_cart))
	ps_cart_in = rbind(ps_cart_in, as.vector(fit_cart))
	index_carts = rbind(index_carts, index_cart)
	rocs_cart = c(rocs_cart, list(roc_cart))
	threshes_cart = c(threshes_cart, thresh_cart)


	##RF
	dt = data.frame(X_train, y_train)
	prior = c(1-one_prob, one_prob)
	names(prior) = c("0", "1")

	model_rf = randomForest(y_train~., data=dt, ntree=500, cutoff=prior)
	fit_rf = predict(object=model_rf, newdata=data.frame(X_train), type="prob")[,2]
	p_rf = predict(object=model_rf, newdata=data.frame(X_test), type="prob")[,2]
	y_rf = predict(object=model_rf, newdata=data.frame(X_test), type="response")

	index_rf = eva(y_test, y_rf, p_rf)
	roc_rf = roc(p_rf, y_test)

	ps_rf_out = rbind(ps_rf_out, as.vector(p_rf))
	ps_rf_in = rbind(ps_rf_in, as.vector(fit_rf))
	index_rfs = rbind(index_rfs, index_rf)
	rocs_rf = c(rocs_rf, list(roc_rf))


	##Logistic(l1)
	fold = 5
	cv = cv.glmnet(x=data.matrix(X_train), y=y_train, family="binomial", type.measure="auc", nfolds=fold)

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


	model = glmnet(x=data.matrix(X_train), y=y_train, family="binomial", alpha=1)
	beta = as.numeric(coef(model, s=lambda)[-1])

	p_la = predict(model, newx=data.matrix(X_test), type="response", s=lambda)
	fit_la = predict(model, newx=data.matrix(X_train), type="response", s=lambda)
	roc_la = roc(p_la, y_test)
	thresh_la = one_prob
	y_la = as.factor(ifelse(p_la>thresh_la, 1, 0))
	index_la = eva(y_test, y_la, p_la)

	ps_la_out = rbind(ps_la_out, as.vector(p_la))
	ps_la_in = rbind(ps_la_in, as.vector(fit_la))
	index_las = rbind(index_las, index_la)
	rocs_la = c(rocs_la, list(roc_la))
	threshes_la = c(threshes_la, thresh_la)


	##Logistic(refit)
	X_temp = X_train[,which(beta!=0)]
	glm_data = data.frame(X_temp, y=y_train)
	glm_model = glm(y~., data=glm_data, family="binomial")

	X_test_temp = data.frame(X_test[,which(beta!=0)])
	colnames(X_test_temp) = colnames(glm_data)[-ncol(glm_data)]
	p_glm= predict(glm_model, newdata=X_test_temp, type="response")
	fit_glm= predict(glm_model, newdata=X_temp, type="response")
	y_glm = as.factor(ifelse(p_glm>one_prob, 1, 0))
	index_glm = eva(y_test, y_glm, p_glm)
	roc_glm = roc(p_glm, y_test)
	thresh_glm = one_prob

	ps_glm_out = rbind(ps_glm_out, as.vector(p_glm))
	ps_glm_in = rbind(ps_glm_in, as.vector(fit_glm))
	index_glms = rbind(index_glms, index_glm)
	rocs_glm = c(rocs_glm, list(roc_glm))
	threshes_glm = c(threshes_glm, thresh_glm)


	###Logistic(only)
	glm_data = data.frame(X_train, y=y_train)
	glm_model = glm(y~., data=glm_data, family="binomial")

	p_glm_only= predict(glm_model, newdata=data.frame(X_test), type="response")
	fit_glm_only= predict(glm_model, newdata=data.frame(X_train), type="response")
	roc_glm_only = roc(p_glm_only, y_test)
	thresh_glm_only = one_prob
	y_glm_only = as.factor(ifelse(p_glm_only>thresh_glm_only, 1, 0))
	index_glm_only = eva(y_test, y_glm_only, p_glm_only)

	ps_glm_only_out = rbind(ps_glm_only_out, as.vector(p_glm_only))
	ps_glm_only_in = rbind(ps_glm_only_in, as.vector(fit_glm_only))
	index_glms_only = rbind(index_glms_only, index_glm_only)
	rocs_glm_only = c(rocs_glm_only, list(roc_glm_only))
	threshes_glm_only = c(threshes_glm_only, thresh_glm_only)


	# save the results in the current loop
	ys_out = rbind(ys_out, y_test)
	ys_in = rbind(ys_in, y_train)
	save(list=savelist, file=file_name)

}

end = proc.time()
end - start



