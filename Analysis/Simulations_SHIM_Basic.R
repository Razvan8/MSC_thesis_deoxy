
libs_path<-file.path("..","libs")
source(file.path(libs_path,'Shim3way_class.R'))
source(file.path(libs_path,'hierarchy_tests.R'))
source(file.path(libs_path,'recover_parameters.R'))



data<- create_basic_dataset()
X<- data$X
y<- data$y$obs

beta_true<- data$beta[-1,]
l1=8
l2=8
l3=4
#print(beta_true)

range_main<-c(1: (l1+l2+l3) )
range_theta<-c( (l1+l2+l3+1) : (l1+l2+l3+l1*(l2+l3)+l2*l3) )
range_psi<-c(  (l1+l2+l3+ 1+ l1*(l2+l3)+l2*l3): (l1+l2+l3+ l1*(l2+l3)+l2*l3+l1*l2*l3) )

beta_main<-beta_true[1:(l1+l2+l3)]
beta_2way<-beta_true[unlist(get_ranges(l1,l2,l3)[2])]
beta_3way<-beta_true[unlist(get_ranges(l1,l2,l3)[3])]
beta_2way_matrix<-get_theta_from_theta_vec_2way3(beta_2way,l1=l1,l2=l2, l3=l3)
beta_3way_table<-get_psi_from_psi_vec3(beta_3way,l1=l1,l2=l2, l3=l3)



beta_main_recovered<- get_all_beta(beta_main, l1=l1, l2=l2, l3=l3, threshold = 1e-4)
beta_2way_recovered<-  get_theta_vec_2way3(  get_all_theta(beta_2way_matrix, l1=l1, l2=l2, l3=l3, threshold = 1e-4), l1=l1+1, l2=l2+1, l3=l3+1)
beta_3way_recovered<- get_psi_vec3( get_all_psi(beta_3way_table, l1=l1, l2=l2, l3=l3, threshold = 1e-4) , l1=l1+1, l2=l2+1, l3=l3+1)


beta_2way_without_gamma<-get_beta_vec_2way(beta = beta_main, l1=l1, l2=l2, l3=l3, gamma=NULL, only_beta = TRUE )
beta_3way_without_gamma<-get_beta_vec_3way(beta_2way = beta_2way, l1=l1, l2=l2, l3=l3, delta = NULL, only_beta = TRUE)


gamma_true<-beta_2way/beta_2way_without_gamma
gamma_true[is.nan(gamma_true)]<-0
delta_true<-beta_3way/beta_3way_without_gamma
delta_true[is.nan(delta_true)]<-0

y_centered<-y-mean(y)


####START LASSO



cv_fit <- cv.glmnet(X, y_centered, alpha = 1)
best_lambda <- cv_fit$lambda.min
cat("best lambda: ",best_lambda)
best_lambda<-0.27
lasso_model <- glmnet(X, y_centered, alpha = 1, intercept = FALSE, standardize = FALSE, lambda=best_lambda)

coefs_lasso<-coefficients(lasso_model)[-1]
beta_main_lasso<-coefs_lasso[range_main]
beta_2way_lasso<-coefs_lasso[range_theta]
beta_3way_lasso<-coefs_lasso[range_psi]
beta_2way_lasso_without_gamma<-get_beta_vec_2way(beta_main_lasso,l1=l1,l2=l2,l3=l3,only_beta = TRUE)
beta_3way_lasso_without_delta<- get_beta_vec_3way(beta_2way_lasso, l1=l1, l2=l2, l3=l3, only_beta = TRUE)

predict_lasso<-predict(lasso_model, s=best_lambda, newx = X)
print(r2(y_centered, predict_lasso))

length(c(beta_main, beta_2way, beta_3way))



######################### RESULTS LASSO #########################################################
print(" results lasso without recovered")
all_beta_functions(beta_main, beta_main_lasso)
all_beta_functions(beta_2way, beta_2way_lasso)
all_beta_functions(beta_3way, beta_3way_lasso)





##### RESULTS LASSO ON RECOVERED PARAMS ##########
beta_2way_lasso_matrix<-get_theta_from_theta_vec_2way3(beta_2way_lasso,l1=l1,l2=l2, l3=l3)
beta_3way_lasso_table<-get_psi_from_psi_vec3(beta_3way_lasso,l1=l1,l2=l2, l3=l3)

##hierarchy tests

test_hierarchy_layer12(beta_main_lasso,beta_2way_lasso_matrix, strong = TRUE)
test_hierarchy_layer23(beta_2way_lasso_matrix, beta_3way_lasso_table, strong = TRUE)

beta_main_lasso_recovered<- get_all_beta(beta_main_lasso, l1=l1, l2=l2, l3=l3, threshold = 1e-4)
beta_2way_lasso_recovered<-  get_theta_vec_2way3(  get_all_theta(beta_2way_lasso_matrix, l1=l1, l2=l2, l3=l3, threshold = 1e-4), l1=l1+1, l2=l2+1, l3=l3+1)
beta_3way_lasso_recovered<- get_psi_vec3( get_all_psi(beta_3way_lasso_table, l1=l1, l2=l2, l3=l3, threshold = 1e-4) , l1=l1+1, l2=l2+1, l3=l3+1)


print("results lasso recovered")
all_beta_functions(beta_main_recovered, beta_main_lasso_recovered)
all_beta_functions(beta_2way_recovered, beta_2way_lasso_recovered)
all_beta_functions(beta_3way_recovered, beta_3way_lasso_recovered)



##PREPARE FOR SHIM ###################################################################################################

beta_hat<-beta_main_lasso
gamma_hat<- beta_2way_lasso/beta_2way_lasso_without_gamma
gamma_hat[is.nan(gamma_hat)]<-0
gamma_hat[!is.finite(gamma_hat)]<-0 #this is 0 in shim case
delta_hat<- beta_3way_lasso/beta_3way_lasso_without_delta
delta_hat[!is.finite(delta_hat)]<-0
delta_hat[is.nan(delta_hat)]<-0

##USE SHIM MODEL #########

lambda_beta<-300
lambda_gamma<-2200
lambda_delta<-700


my_shim<-SHIM_3way(X=X, y=y, beta_init = beta_hat, gamma_init = gamma_hat, delta_init = delta_hat, l1=l1, l2=l2, l3=l3, scale = FALSE)
fitted<-my_shim$fit(X=X, y=y, lambda_beta = lambda_beta, 
                    lambda_gamma = lambda_gamma, lambda_delta = lambda_delta, w_beta = 1, w_gamma = 1, w_delta = 1, tol=1e-2)

my_shim$R2_score(self=fitted, X_new=X, y_true=y )


beta_all_shim<-fitted$beta_all
beta_main_shim<-beta_all_shim[range_main]
beta_2way_shim<-beta_all_shim[range_theta]
beta_3way_shim<-beta_all_shim[range_psi]



all_beta_functions(beta_main, beta_main_shim)
all_beta_functions(beta_2way, beta_2way_shim)
all_beta_functions(beta_3way, beta_3way_shim)

##hierarchy tests
beta_2way_shim_matrix<-get_theta_from_theta_vec_2way3(beta_2way_shim,l1=l1,l2=l2, l3=l3)
beta_3way_shim_table<-get_psi_from_psi_vec3(beta_3way_shim,l1=l1,l2=l2, l3=l3)

test_hierarchy_layer12(beta_main_shim,beta_2way_shim_matrix, strong = TRUE)
test_hierarchy_layer23(beta_2way_shim_matrix, beta_3way_shim_table, strong = TRUE)



### SHIM FOR RECOVERED DATA #############
beta_2way_shim_matrix<-get_theta_from_theta_vec_2way3(beta_2way_shim,l1=l1,l2=l2, l3=l3)
beta_3way_shim_table<-get_psi_from_psi_vec3(beta_3way_shim,l1=l1, l2=l2, l3=l3)

beta_main_shim_recovered<- get_all_beta(beta_main_shim, l1=l1, l2=l2, l3=l3, threshold = 1e-4)
beta_2way_shim_recovered<-  get_theta_vec_2way3( get_all_theta(beta_2way_shim_matrix, l1=l1, l2=l2, l3=l3, threshold = 1e-4), l1=l1+1, l2=l2+1, l3=l3+1)
beta_3way_shim_recovered<- get_psi_vec3( get_all_psi(beta_3way_shim_table, l1=l1, l2=l2, l3=l3, threshold = 1e-4) , l1=l1+1, l2=l2+1, l3=l3+1)


print("results shim recovered")
all_beta_functions(beta_main_recovered, beta_main_shim_recovered)
all_beta_functions(beta_2way_recovered, beta_2way_shim_recovered)
all_beta_functions(beta_3way_recovered, beta_3way_shim_recovered)
















####### SAME ANALYSIS BUT DIFFERENT LAMBDA (smaller - get aroung 80% of zeros) ###################
####START LASSO



cv_fit <- cv.glmnet(X, y_centered, alpha = 1)

best_lambda<-0.08
lasso_model <- glmnet(X, y_centered, alpha = 1, intercept = FALSE, standardize = FALSE, lambda=best_lambda)

coefs_lasso<-coefficients(lasso_model)[-1]
beta_main_lasso<-coefs_lasso[range_main]
beta_2way_lasso<-coefs_lasso[range_theta]
beta_3way_lasso<-coefs_lasso[range_psi]
beta_2way_lasso_without_gamma<-get_beta_vec_2way(beta_main_lasso,l1=l1,l2=l2,l3=l3,only_beta = TRUE)
beta_3way_lasso_without_delta<- get_beta_vec_3way(beta_2way_lasso, l1=l1, l2=l2, l3=l3, only_beta = TRUE)

predict_lasso<-predict(lasso_model, s=best_lambda, newx = X)
print(r2(y_centered, predict_lasso))

length(c(beta_main, beta_2way, beta_3way))

######################### RESULTS LASSO #########################################################
print(" results lasso without recovered")
all_beta_functions(beta_main, beta_main_lasso)
all_beta_functions(beta_2way, beta_2way_lasso)
all_beta_functions(beta_3way, beta_3way_lasso)


##hierarchy tests

beta_2way_lasso_matrix<-get_theta_from_theta_vec_2way3(beta_2way_lasso,l1=l1,l2=l2, l3=l3)
beta_3way_lasso_table<-get_psi_from_psi_vec3(beta_3way_lasso,l1=l1,l2=l2, l3=l3)

test_hierarchy_layer12(beta_main_lasso,beta_2way_lasso_matrix, strong = TRUE)
test_hierarchy_layer23(beta_2way_lasso_matrix, beta_3way_lasso_table, strong = TRUE)


##### RESULTS LASSO ON RECOVERED PARAMS ##########

beta_main_lasso_recovered<- get_all_beta(beta_main_lasso, l1=l1, l2=l2, l3=l3, threshold = 1e-6)
beta_2way_lasso_recovered<-  get_theta_vec_2way3(  get_all_theta(beta_2way_lasso_matrix, l1=l1, l2=l2, l3=l3, threshold = 1e-6), l1=l1+1, l2=l2+1, l3=l3+1)
beta_3way_lasso_recovered<- get_psi_vec3( get_all_psi(beta_3way_lasso_table, l1=l1, l2=l2, l3=l3, threshold = 1e-6) , l1=l1+1, l2=l2+1, l3=l3+1)


#beta_main_lasso
#beta_main_lasso_recovered

print("results lasso recovered")
all_beta_functions(beta_main_recovered, beta_main_lasso_recovered)
all_beta_functions(beta_2way_recovered, beta_2way_lasso_recovered)
all_beta_functions(beta_3way_recovered, beta_3way_lasso_recovered)


test_hierarchy_layer12(beta_main_lasso_recovered, get_theta_from_theta_vec_2way3( beta_2way_lasso_recovered, l1=l1+1, l2=l2+1, l3=l3+1 ), strong =TRUE)
test_hierarchy_layer23( get_theta_from_theta_vec_2way3( beta_2way_lasso_recovered, l1=l1+1, l2=l2+1, l3=l3+1), 
                        get_psi_from_psi_vec3(beta_3way_lasso_recovered,l1=l1+1,l2=l2+1, l3=l3+1), strong =TRUE)



##PREPARE FOR SHIM ###################################################################################################

beta_hat<-beta_main_lasso
gamma_hat<- beta_2way_lasso/beta_2way_lasso_without_gamma
gamma_hat[is.nan(gamma_hat)]<-0
gamma_hat[!is.finite(gamma_hat)]<-0 #this is 0 in shim case
delta_hat<- beta_3way_lasso/beta_3way_lasso_without_delta
delta_hat[!is.finite(delta_hat)]<-0
delta_hat[is.nan(delta_hat)]<-0

##USE SHIM MODEL #########

lambda_beta<-80
lambda_gamma<-1000
lambda_delta<-1000


my_shim<-SHIM_3way(X=X, y=y, beta_init = beta_hat, gamma_init = gamma_hat, delta_init = delta_hat, l1=l1, l2=l2, l3=l3, scale = FALSE)
fitted<-my_shim$fit(X=X, y=y, lambda_beta = lambda_beta, 
                    lambda_gamma = lambda_gamma, lambda_delta = lambda_delta, w_beta = 1, w_gamma = 1, w_delta = 1, tol=8e-4)

my_shim$R2_score(self=fitted, X_new=X, y_true=y )


beta_all_shim<-fitted$beta_all
beta_main_shim<-beta_all_shim[range_main]
beta_2way_shim<-beta_all_shim[range_theta]
beta_3way_shim<-beta_all_shim[range_psi]



all_beta_functions(beta_main, beta_main_shim)
all_beta_functions(beta_2way, beta_2way_shim)
all_beta_functions(beta_3way, beta_3way_shim)

##hierarchy tests
beta_2way_shim_matrix<-get_theta_from_theta_vec_2way3(beta_2way_shim,l1=l1,l2=l2, l3=l3)
beta_3way_shim_table<-get_psi_from_psi_vec3(beta_3way_shim,l1=l1,l2=l2, l3=l3)

test_hierarchy_layer12(beta_main_shim,beta_2way_shim_matrix, strong = TRUE)
test_hierarchy_layer23(beta_2way_shim_matrix, beta_3way_shim_table, strong = TRUE)



### SHIM FOR RECOVERED DATA #############
beta_2way_shim_matrix<-get_theta_from_theta_vec_2way3(beta_2way_shim,l1=l1,l2=l2, l3=l3)
beta_3way_shim_table<-get_psi_from_psi_vec3(beta_3way_shim,l1=l1, l2=l2, l3=l3)

beta_main_shim_recovered<- get_all_beta(beta_main_shim, l1=l1, l2=l2, l3=l3, threshold = 1e-4)
beta_2way_shim_recovered<-  get_theta_vec_2way3( get_all_theta(beta_2way_shim_matrix, l1=l1, l2=l2, l3=l3, threshold = 1e-4), l1=l1+1, l2=l2+1, l3=l3+1)
beta_3way_shim_recovered<- get_psi_vec3( get_all_psi(beta_3way_shim_table, l1=l1, l2=l2, l3=l3, threshold = 1e-4) , l1=l1+1, l2=l2+1, l3=l3+1)


print("results shim recovered")
all_beta_functions(beta_main_recovered, beta_main_shim_recovered)
all_beta_functions(beta_2way_recovered, beta_2way_shim_recovered)
all_beta_functions(beta_3way_recovered, beta_3way_shim_recovered)

test_hierarchy_layer12(beta_main_shim_recovered, get_theta_from_theta_vec_2way3( beta_2way_shim_recovered, l1=l1+1, l2=l2+1, l3=l3+1 ), strong =TRUE)
test_hierarchy_layer23( get_theta_from_theta_vec_2way3( beta_2way_shim_recovered, l1=l1+1, l2=l2+1, l3=l3+1), 
                        get_psi_from_psi_vec3(beta_3way_shim_recovered,l1=l1+1,l2=l2+1, l3=l3+1), strong =TRUE)




