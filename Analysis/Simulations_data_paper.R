####Sparse hierarchical simulations####

library(glmnet)

libs_path<-file.path("..","libs")
source(file.path(libs_path,'Create_synthetic_datasets.R'))
source(file.path(libs_path,'WeakHierNet_Class_corrected_unscaled.R'))
source(file.path(libs_path,'WeakHierNetSeq23_3way.R'))



data=create_hier_dataset_paper_3way()
X<-data$X
y<-data$y$obs
beta_trues<-data$beta[-1,]  ##without intercept






l1=9
l2=9
l3=4



#y
#beta_trues

range_main<-c(1: (l1+l2+l3) )
range_theta<-c( (l1+l2+l3+1) : (l1+l2+l3+l1*(l2+l3)+l2*l3) )
range_psi<-c(  (l1+l2+l3+ 1+ l1*(l2+l3)+l2*l3): (l1+l2+l3+ l1*(l2+l3)+l2*l3+l1*l2*l3) )

colnames(X)[range_main]

beta_true<-beta_trues[range_main]
theta_true<-beta_trues[range_theta]
psi_true<-beta_trues[range_psi]

#theta_true

print(dim(X))
print(dim(y))


dim(X)


X_only_main<-X[,range_main]

X_2way<-X[,range_theta]
X_3way<-X[,range_psi]
y_all<-as.vector(y[,1])






################################### ANALYSIS WITH LASSO ##############################################################

####Lasso only on main effects
#lasso_model <- glmnet(X[,range_main], y_all, alpha = 1, lambda=0.3)
#coefs<-coefficients(lasso_model)[-1]
#coefs
#sum(coefs==0)

#all_beta_functions(beta_true, coefs[range_main])

#print(beta_true)


sum(beta_trues==0)

lasso_model <- glmnet(X, y_all, alpha = 1, lambda=0.25)
coefs<-coefficients(lasso_model)[-1]
#coefs
#sum(coefs==0)

print(coefs[range_theta])

all_beta_functions(beta_true, coefs[range_main])
all_beta_functions(theta_true, coefs[range_theta])
all_beta_functions(psi_true, coefs[range_psi])

sum(coefs==0)

#####

beta_init_lasso<-coefs[range_main]
theta_init<-get_theta_hat_from_vec3(coefs[range_theta],l1=l1,l2=l2,l3=l3)

print(theta_init)
#theta_init<-matrix(0, ncol = dim(X_only_main)[2]



#######ANALYSIS SEQ1-2 + SEQ 2-3 #######################




#####

beta_init_lasso<-coefs[range_main]
theta_init<-get_theta_hat_from_vec3(coefs[range_theta],l1=l1,l2=l2,l3=l3)
coeffs_main<-coefs[range_main]

#theta_init<-matrix(0, nrow=dim(X_only_main)[2], ncol = dim(X_only_main)[2])




#dim(X_only_main)[2]
#l1+l2+l3


########################## ANALYSIS WITH SEQ1-2 and SEQ2-3 ###########################################################



#all_beta_functions(theta_true, c(fitted$vec_theta_hat))

#fitted$theta_hat

#print(fitted$vec_theta_hat[0:30])
#print(theta_true[0:30])

#print(coeffs_main)
#print(fitted$theta_hat)

#print(length(range_theta))















###### ANALYSIS WITH my library   #########################################

lmd<-50
t<-1e-5

beta_init_lasso_plus<- beta_init_lasso
beta_init_lasso_plus[beta_init_lasso_plus<0]<-0

beta_init_lasso_minus<- -beta_init_lasso
beta_init_lasso_plus[beta_init_lasso_plus<0]<-0






source(file.path(libs_path,'WeakHierNet_Class_corrected_unscaled.R'))
myWeakHierNet<-WeakHierNetUnscaled (X=X_only_main, Beta_plus_init=beta_init_lasso_plus, Beta_minus_init=beta_init_lasso_minus, 
                                    Theta_init=theta_init, y=y_all, lambda=lmd, t=t, tol=1e-8, 
                                    max_iter=10000, eps=1e-8, center = FALSE, standard_form=TRUE)  #Increase max iter if needed or decrease tol 

# Fit the model
fitted=myWeakHierNet$fit(X=X_only_main,Beta_plus_init=beta_init_lasso_plus, Beta_minus_init = beta_init_lasso_minus, 
                         Theta_init=theta_init, y=y_all, lambda=lmd, t=t, tol=1e-8, 
                         max_iter=10000, eps=1e-8,l1=l1,l2=l2,l3=l3 )

print(fitted$Beta_hat_plus-fitted$Beta_hat_minus)

print(fitted$Beta_hat_plus)
print(fitted$Beta_hat_minus)

print(fitted$Beta_hat_plus-fitted$Beta_hat_minus)
#print(fitted$vec_theta_hat)
print(fitted$Theta_hat)
print(get_vec_theta_hat3(fit$th,l1=l1,l2=l2,l3=l3))

#sum(fitted$vec_theta_hat==0)
#sum(fitted$Beta_hat_plus-fitted$Beta_hat_minus==0)
#print(myWeakHierNet$R2_score(fitted,X_only_main,y_all))

print(fitted$vec_theta_hat)
print(beta_trues[range_theta])
print(get_vec_theta_hat3(fit$th,l1=l1,l2=l2,l3=l3))

print(fit$th)

###RESULTS FOR BETA AND THETA#####
print("My library")
all_beta_functions(beta_true, fitted$Beta_hat_plus-fitted$Beta_hat_minus)
all_beta_functions(theta_true, fitted$vec_theta_hat)

sum(fitted$vec_theta_hat==0)+sum(fitted$Beta_hat_plus-fitted$Beta_hat_minus ==0)
sum(beta_init_lasso_plus-beta_init_lasso_minus ==0) +sum(coefs[range_theta]==0)
coefs[range_theta]==get_vec_theta_hat3(theta_init,l1=l1,l2=l2,l3=l3)

print(fitted$Theta_hat)
fitted$vec_theta_hat
beta_trues[range_theta]

###HIERNET LIBRARY##############################################
#print("-----Hiernet library-----")

fit=hierNet(X_only_main, y_all, lam=30, diagonal = FALSE, stand.main=FALSE,tol=1e-8)
#fit$th
predicted=predict(fit,X_only_main)
print(r2(y_all, predicted))

#fit$bp-fit$bn
#print(sum( get_vec_theta_hat3(fit$th,l1=l1, l2=l2, l3=l3)==0))
print("hiernet")
all_beta_functions(beta_true, c(fit$bp-fit$bn) )
all_beta_functions(theta_true, get_vec_theta_hat3(fit$th,l1=l1, l2=l2, l3=l3))
#get_vec_theta_hat3(fit$th,l1=l1, l2=l2, l3=l3)[1:40]
#theta_true[1:40]
fit$th





#######

#Residuals
predicted<-myWeakHierNet$predict(fitted, X_only_main)
r_2way<-y_all-predicted


##Fit in on residuals

lasso_model <- glmnet(X_3way, r_2way, alpha = 1, lambda=0.05)
coefficients(lasso_model)

psi_init<-get_psi_from_psi_vec3(as.vector(coefficients(lasso_model)[-1]),l1=l1,l2=l2,l3=l3) ## psi, get intercept out
print(dim(psi_init))
print(dim(X_3way))
theta_bound<-fitted$Theta_hat +t(fitted$Theta_hat)*1 ##bound
lambda<-70

source(file.path(libs_path,'WeakHierNetSeq23_3way.R'))
t<-1e-1
myWeakHierNet_seq3 <- WeakHierNet_seq3(X=X_3way, psi_init=psi_init, y=r_2way, theta_bound=theta_bound, lambda=lambda, t=t, tol=1e-9, max_iter=5000, eps=1e-8,
                                    l1=l1,l2=l2,l3=l3, scale = FALSE)

# Fit the model
fitted3=myWeakHierNet_seq3$fit(X=X_3way, psi_init=psi_init, y=r_2way, theta_bound=theta_bound, lambda=lambda, t=t, tol=1e-9, max_iter=5000, 
                            eps=1e-8,l1=l1,l2=l2,l3=l3)

r2(r_2way,myWeakHierNet_seq3$predict(fitted3,X_3way))
#print(fitted$vec_psi_hat)

sum(fitted3$vec_psi_hat==0)
dim(fitted3$vec_psi_hat)

fitted3$psi_hat[8,15,22]

###RESULTS FOR PSI#####

all_beta_functions(psi_true, c(fitted3$vec_psi_hat))
vec_nulled<-fitted3$vec_psi_hat
vec_nulled[abs(vec_nulled)<1e-1]<-0
all_beta_functions(psi_true, vec_nulled)
all_beta_functions(psi_true,coefs[range_psi])
vec_nulled
#######

fitted$vec_psi_hat



########################## ANALYSIS WITH SEQ1-2 and SEQ2-3 ###########################################################


