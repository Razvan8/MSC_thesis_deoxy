####Sparse hierarchical simulations####

library(glmnet)

libs_path<-file.path("..","libs")
source(file.path(libs_path,'Create_synthetic_datasets.R'))
source(file.path(libs_path,'WeakHierNet_Class_corrected_unscaled.R'))
source(file.path(libs_path,'WeakHierNetSeq23_3way.R'))



data=create_sparse_hier_dataset()
X<-data$X
y<-data$y
print(dim(X))
print(dim(y))

colnames(X)[92]

X_only_main<-X[,c(1:16)]
X_3way<-X[,c(92:199)]
y_all<-as.vector(y[,1])
y_all


###HIERNET LIBRARY
print("-----Hiernet library-----")
##On train
fit=hierNet(X_only_main, y_all, lam=10, diagonal = FALSE, stand.main=FALSE,tol=1e-6)
#fit
#y_pred_train=predict(fit,X_only_main)
#sum(fit$bp+fit$bn)
#$sum(abs(fit$th))
print(get_vec_theta_hat3(fit$th,l1=4, l2=9, l3=3))
#print(fit$bp-fit$bn)



#####
####TRUTH
#beta[1]=1 #A1
#beta[4]=4 #A4
#beta[9]=9 #B5
#beta[13]=-13 #B9
#beta[16]=-1 #C3

#2 way
#beta[21]=2 #A1:B5
#beta[25]=-2 #A1:B9
#beta[91]=9 #B9:C3

#3way
#beta[106]=2  #A1B5C3  sau 106-91=15
#beta[105]=-1  #A1B5C2             =14 
#beta[128]=-2  #A1B9C3  sau 128-91=37



#####
lmd<-10
t<-1e-3

myWeakHierNet<-WeakHierNetUnscaled (X=X_only_main, Beta_plus_init=matrix(0,nrow = dim(X_only_main)[2], ncol = 1), Beta_minus_init=matrix(0,nrow = dim(X_only_main)[2], ncol = 1), 
                                    Theta_init=matrix(0, ncol = dim(X_only_main)[2], nrow = dim(X_only_main)[2]), y=y_all, lambda=lmd, t=t, tol=1e-8, 
                                    max_iter=10000, eps=1e-8)  #Increase max iter if needed or decrease tol 

# Fit the model
fitted=myWeakHierNet$fit(X=X_only_main, Beta_plus_init=matrix(0,nrow = dim(X_only_main)[2], ncol = 1), Beta_minus_init=matrix(0,nrow = dim(X_only_main)[2], ncol = 1), 
                         Theta_init=matrix(0, ncol = dim(X_only_main)[2], nrow = dim(X_only_main)[2]), y=y_all, lambda=lmd, t=t, tol=1e-8, 
                         max_iter=10000, eps=1e-8,l1=4,l2=9,l3=3 )

print(fitted$Beta_hat_plus-fitted$Beta_hat_minus)
print(fitted$vec_theta_hat)
print(fitted$Theta_hat)
print(myWeakHierNet$R2_score(fitted,X_only_main,y_all))

#Residuals
predicted<-myWeakHierNet$predict(fitted, X_only_main)
r_2way<-y_all-predicted


##Fit in on residuals

lasso_model <- glmnet(X_3way, r_2way, alpha = 1, lambda=1)
coefficients(lasso_model)

psi_init<-get_psi_from_psi_vec3(as.vector(coefficients(lasso_model)[-1]),l1=4,l2=9,l3=3) ## psi, get intercept out
theta_bound<-fitted$Theta_hat +t(fitted$Theta_hat) ##bound
lambda<-10


myWeakHierNet_seq3 <- WeakHierNet_seq3(X=X_3way, psi_init=psi_init, y=r_2way, theta_bound=theta_bound, lambda=lambda, t=t, tol=1e-8, max_iter=5000, eps=1e-8,
                                    l1=4,l2=9,l3=3, scale = FALSE)

# Fit the model
fitted=myWeakHierNet_seq3$fit(X=X_3way, psi_init=psi_init, y=r_2way, theta_bound=theta_bound, lambda=lambda, t=t, tol=1e-8, max_iter=5000, 
                            eps=1e-8,l1=4,l2=9,l3=3)

r2(r_2way,myWeakHierNet_seq3$predict(fitted,X_3way))


print(fitted$vec_psi_hat)







