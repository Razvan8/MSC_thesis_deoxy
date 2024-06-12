## FUNCTIONS USED FOR THE CLASS

assert <- function(condition, message) {
  if (!condition) stop(message)
}




get_xx.all<-function(X,l1,l2,l3)
{xx.all<- matrix(0, nrow=dim(X)[1], ncol=l1+l2+l3+l1*(l2+l3) +l2*l3 )
xx.all[,c(1:(l1+l2+l3))]<-X
counter<-l1+l2+l3+1

for (i in c(1:l1)){ #ab ac
  for (j in c ( (l1+1): (l1+l2+l3) ) ){
    print(counter)
    xx.all[,counter]<-X[,i]*X[,j]
    counter<-counter+1}}

for (i in c((l1+1): (l1+l2))){ #bc
  for (j in  c( (l1+l2+1): (l1+l2+l3) ) ){
    xx.all[,counter]<-X[,i]*X[,j]
    counter<-counter+1
    }}

assert(counter== l1+l2+l3+ l1*(l2+l3)+l2*l3+1)
return(xx.all)
}

get_xxx.all<-function(X,l1,l2,l3)
{xxx.all<- matrix(0, nrow=dim(X)[1], ncol=l1+l2+l3+l1*(l2+l3) +l2*l3 + l1*l2*l3)
xxx.all[,c(1:(l1+l2+l3+ l1*(l2+l3)+l2*l3))]<-get_xx.all(X=X, l1=l1, l2=l2, l3=l3)
counter<-l1+l2+l3+l1*(l2+l3)+l2*l3+1

for (i in c(1:l1)){ #abc
  for (j in c ( (l1+1): (l1+l2) ) ){
    for (k in c ( (l1+l2+1): (l1+l2+l3) ) ){
      
      xxx.all[,counter]<-X[,i]*X[,j]*X[,k]
      counter<-counter+1}}}
assert(counter==l1+l2+l3+l1*(l2+l3)+l2*l3+l1*l2*l3+1)
return(xxx.all)
}



get_beta_vec_2way<-function(beta,l1,l2,l3)
{beta_vec2way<- array(0, dim = l1*(l2+l3) +l2*l3 )
counter<-1
  
  for (i in c(1:l1)){ #ab ac
  for (j in c ( (l1+1): (l1+l2+l3) ) ){
    beta_vec2way[counter]<-beta[i]*beta[j]
    counter<-counter+1}}

  for (i in c((l1+1): (l1+l2))){ #bc
    for (j in  c( (l1+l2+1): (l1+l2+l3) ) ){
      beta_vec2way[counter]<-beta[i]*beta[j]
      counter<-counter+1}}

assert(counter==l1*(l2+l3)+l2*l3+1)
return(beta_vec2way)
}


get_beta_vec_3way<-function(beta,l1,l2,l3)
{beta_vec3way<- array(0, dim = l1*l2*l3 )
counter<-1

for (i in c(1:l1)){ #ab ac
  for (j in c ( (l1+1): (l1+l2) ) ){
    for (k in c ( (l1+l2+1): (l1+l2+l3) ) ){
    
    beta_vec3way[counter]<-beta[i]*beta[j]*beta[k]
    counter<-counter+1}}}


assert(counter==l1*l2*l3+1)

return(beta_vec3way)

}





get_ranges<-function(l1,l2,l3)
{range_main<-c(1: (l1+l2+l3) )
range_theta<-c( (l1+l2+l3+1) : (l1+l2+l3+l1*(l2+l3)+l2*l3) )
range_psi<-c(  (l1+l2+l3+ 1+ l1*(l2+l3)+l2*l3): (l1+l2+l3+ l1*(l2+l3)+l2*l3+l1*l2*l3) )
return(list(range_main, range_theta, range_psi))}
  


mains_contribution<-function(X, beta, l1,l2,l3)
{ range_main<-unlist(get_ranges(l1,l2,l3)[1])
  mains_contrib<-X[,range_main]%*%beta
  return(mains_contrib)}

two_ways_contribution<-function(X, gamma_vec, beta_vec_2way,l1,l2,l3) ##assumes gamma_vec is in same order as X
{range_2ways<-unlist(get_ranges(l1,l2,l3)[2])
#print(dim())
 two_ways_contrib<- X[,range_2ways]%*%(beta_vec_2way*gamma_vec) ##last multiplication should be elementwise
 return(two_ways_contrib)}

three_ways_contribution<-function(X, delta_vec, beta_vec_3way, l1,l2,l3) ##assumes gamma_vec is in same order as X
{range_3ways<-unlist(get_ranges(l1,l2,l3)[3])
three_ways_contrib<-X[,range_3ways]%*%(beta_vec_3way*delta_vec) ##last multiplication should be elementwise
return(three_ways_contrib)}



### COMPUTE Q

compute_Q<-function()


###RELATIVE DIFFERENCE

compute_relative_dif<-function(Q_old, Q_new)
{rel_dif<- abs(Q_old-Q_new)/abs(Q_old)
return(rel_dif)}



##### UPDATE DELTA FUNCTION #####

update_delta<-function(beta_hat, gamma_hat, delta_hat, lambda_delta)
{
  
  
  
}




##### UPDATE GAMMA FUNCTION #####

update_gamma<-function(beta_hat, gamma_hat, delta_hat, lambda_gamma)
{
  
  
  
}



##### UPDATE BETA FUNCTION #####

update_beta<-function(beta_hat, gamma_hat, delta_hat, lambda_beta)
{
  
  
  
}






#



l1=2
l2=1
l3=1

X=matrix(1, nrow=4, ncol=4)
X[2,4]=0
beta=array(c(1,1,1,2))
beta_vec_2way<-get_beta_vec_2way(beta,l1,l2,l3)
beta_vec_3way<-get_beta_vec_3way(beta,l1,l2,l3)
gamma_vec_2way<-array(c(1,2,3,4,0),dim=length(beta_vec_2way))
delta_vec_3way<-array(c(2,3),dim=length(beta_vec_3way))
xx.all<-get_xx.all(X,l1,l2,l3)
xxx.all<-get_xxx.all(X,l1,l2,l3)



mains_contribution(xxx.all,beta,l1,l2,l3)
two_ways_contribution(xxx.all, gamma_vec_2way, beta_vec_2way,l1,l2,l3 )
three_ways_contribution(xxx.all, beta_vec_3way = beta_vec_3way, delta_vec = delta_vec_3way, l1 , l2, l3 )

## SHIM CLASS for 3 way
SHIM_3way<-function(X, beta_init, gamma_init, delta_init,l1=36,l2=3,l3=4, scale=FALSE)
  
{
  self = list()
  
  self$beta_hat <- beta_init #matrix form
  self$gamma_hat<- gamma_init
  self$delta_hat<-delta_init
  
  self$means_X<-colMeans(X)
  self$stds_X<-apply(X,2,sd)
  self$mean_y<-mean(y)
  self$scale=scale
  
  self$l1=l1
  self$l2=l2
  self$l3=l3
  
  
  
  fit<-function(lambda_beta, lambda_gamma, lambda_delta, w_beta=NULL, w_gamma=NULL, w_delta=NULL,  tol=1e-6, max_iter=50)
  {## STEP 0 (STANDARDIZE)
    if (self$scale == TRUE)
    {      print('was scaled')
      X <- scale(X)}#standardize X
    y <- scale(y, center = TRUE, scale = FALSE) #center y # MAYBE ALREADY SCALED ACTUALLY???????????????????????
    
  
    ## STEP 1 (INIT BETA AND GAMMA AND DELTA)
    beta_hat<-self$beta_hat
    gamma_hat<-self$gamma_hat
    delta_hat<-self$delta_hat
    Q_old<-1e100
    
    for (i in c(1:max_iter))
    {    
      ## STEP 2 (UPDATE DELTA)
      delta_hat<- update_delta(beta_hat, gamma_hat, delta_hat, lambda_delta)
      
      ## STEP 3 (UPDATE GAMMA)
      gamma_hat<- update_gamma(beta_hat, gamma_hat, delta_hat, lambda_delta)
      
      ## STEP 4 (UPDATE BETA)
      beta_hat<- update_beta(beta_hat, gamma_hat, delta_hat, lambda_delta)
      
      ## STEP 5 (COMPUTE REL_DIF)
      Q_new<-compute_Q(.....)
      rel_dif<-compute_rel_dif(Q_old=Q_old, Q_new=Q_new)
      
      
      
      ## STEP 6 (CHECK IF SMALL ENOUGH TO STOP)
      if (abs(rel_dif)<=tol){
        self$beta_hat<-beta_hat
        self$gamma_hat<-gamma_hat
        self$delta_hat<-delta_hat
        return (list("beta_hat" = self$beta_hat, "gamma_hat" = self$gamma_hat , "delta_hat" = self$delta_hat)) }
      Q_old<-Q_new #UPDATE Q_old
      }
      cat("It has not converged. The relative difference between last 2 residuals is:", compute_rel_dif(Q_old=Q_old,Q_new=Q_new) )
      self$beta_hat<-beta_hat
      self$gamma_hat<-gamma_hat
      self$delta_hat<-delta_hat
      return (list("beta_hat" = self$beta_hat, "gamma_hat" = self$gamma_hat , "delta_hat" = self$delta_hat)) 
      }

    
    
}
  
  
  

