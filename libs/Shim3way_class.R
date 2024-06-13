## FUNCTIONS USED FOR THE CLASS

assert <- function(condition, message) {
  if (!condition) stop(message)
}

get_ranges<-function(l1,l2,l3)
{range_main<-c(1: (l1+l2+l3) )
range_theta<-c( (l1+l2+l3+1) : (l1+l2+l3+l1*(l2+l3)+l2*l3) )
range_psi<-c(  (l1+l2+l3+ 1+ l1*(l2+l3)+l2*l3): (l1+l2+l3+ l1*(l2+l3)+l2*l3+l1*l2*l3) )
return(list(range_main, range_theta, range_psi))}


get_range3<- function(x,l1=36,l2=3,l3=4)
{if (x<=l1)
{return(c(1:l1))}
  if (x<=l1+l2)
  {return(c( (l1+1) : (l1+l2) ))}
  
  return(c( (l1+l2+1) : (l1+l2+l3) ))
  
}




##position in matrix form to position in vector form 2way
matrix_position_to_vector_index_2way<- function(position_tuple, l1,l2,l3) ## takes into account / works only for possible combinations!!!!
{ x<-position_tuple[1]
  y<-position_tuple[2]
  
  range_x<-get_range3(x,l1=l1,l2=l2,l3=l3)
  range_y<-get_range3(y,l1=l1,l2=l2,l3=l3)
  

  assert(x<= l1+l2, "x should be <=l1+l2")
  assert(x<y, "x<y")
  assert(y>l1, 'y should be >l1')

  
  if( all(range_x == c(1:l1)) ==TRUE ) #ab or ac
  { 
    position_vector<- (x-1)*(l2+l3) +(y-l1)  }
  
  
  if( all ( range_x == c( (l1+1): (l1+l2) ) ) == TRUE )  #bc
  {position_vector<-l1*(l2+l3) + (x-l1-1)*l3 + y- (l1+l2)  } 
  return(position_vector)
  
}


#position in 3 dim table to vector index: works only for possible combinations
table_position_to_vector_index3<- function(position_tuple, l1,l2,l3) ## takes into account / works only for possible combinations!!!!
{
  
  l12<-l1+l2
  
  x<-position_tuple[1]
  y<-position_tuple[2]
  z<-position_tuple[3]
  
  assert(x<=l1, "x should be <=l1")
  assert(y<=l1+l2, "y should be <=l1+l2")
  assert(y>l1, "y should be >l1")
  assert(z>l1+l2, 'z should be >l1+l2')
  
  position_psi<-(x-1)*l2*l3 + (y-l1-1)*l3 + (z-l12) 
  
  return(position_psi)
  
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



get_beta_vec_2way<-function(beta,l1,l2,l3, gamma, only_beta = FALSE)
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
if (only_beta==FALSE)
{beta_vec2way<-beta_vec2way*gamma}
return(beta_vec2way)
}


get_beta_vec_3way<-function(beta_2way,l1,l2,l3, delta, only_beta=FALSE) ##
# beta_2way should be final beta_2way; 
#only_beta means product of beta_2ways (gamma included) without delta
  
{beta_vec3way<- array(0, dim = l1*l2*l3 )
counter<-1

#Iterate over possible positions
for (i in c(1:l1)){ #ab ac
  for (j in c ( (l1+1): (l1+l2) ) ){
    for (k in c ( (l1+l2+1): (l1+l2+l3) ) ){
    
      beta_vec3way[counter]<-beta_2way[matrix_position_to_vector_index_2way(position_tuple = c(i,j),l1=l1,l2=l2,l3=l3)]*
                           beta_2way[matrix_position_to_vector_index_2way(position_tuple = c(i,k),l1=l1,l2=l2,l3=l3)]*
                           beta_2way[matrix_position_to_vector_index_2way(position_tuple = c(j,k),l1=l1,l2=l2,l3=l3)]
                           
    counter<-counter+1}}}

if (only_beta == FALSE)
{beta_vec3way<-beta_vec3way*delta}


assert(counter==l1*l2*l3+1)

return(beta_vec3way)

}




mains_contribution<-function(X, beta_main, l1,l2,l3)
{ range_main<-unlist(get_ranges(l1,l2,l3)[1])
  mains_contrib<-X[,range_main]%*%beta_main
  return(mains_contrib)}

two_ways_contribution<-function(X, gamma_vec, beta_vec_2way,l1,l2,l3) ##assumes gamma_vec is in same order as X, beta vec is only beta
{range_2ways<-unlist(get_ranges(l1,l2,l3)[2])
#print(dim())
 two_ways_contrib<- X[,range_2ways]%*%(beta_vec_2way*gamma_vec) ##last multiplication should be elementwise
 return(two_ways_contrib)}

three_ways_contribution<-function(X, delta_vec, beta_vec_3way, l1,l2,l3) ##assumes gamma_vec is in same order as X and beta_vec_3way only prod of beta2way
{range_3ways<-unlist(get_ranges(l1,l2,l3)[3])
three_ways_contrib<-X[,range_3ways]%*%(beta_vec_3way*delta_vec) ##last multiplication should be elementwise
return(three_ways_contrib)}



### g function

g_normal<-function(X, beta, gamma_vec, delta_vec,l1,l2,l3, already_multiplied=TRUE) #bet_2way is withoug gamma, beta_3way withour delta only
  #already multiplied=True means beta already has gamma delta inside
{beta_main<-beta[unlist(get_ranges(l1,l2,l3)[1])]
 beta_2way<-beta[unlist(get_ranges(l1,l2,l3)[2])]
 beta_3way<-beta[unlist(get_ranges(l1,l2,l3)[3])]
 if (already_multiplied==TRUE)
 {gamma_vec<-array(1, dim=length(gamma_vec))
  delta_vec<-array(1, dim=length(delta_vec))}
 result<-mains_contribution(X=X, beta_main = beta_main, l1=l1, l2=l2, l3=l3)+ 
         two_ways_contribution(X=X, gamma_vec=gamma_vec, beta_vec_2way=beta_2way,l1=l1, l2=l2, l3=l3)+
         three_ways_contribution(X=X, delta_vec = delta_vec, beta_vec_3way = beta_3way,l1=l1, l2=l2, l3=l3)
 cat("g:", result)
 return(result)
}



##penalty for 1 vector
get_penalty<-function(vector, weights, lambda, already_weighted=TRUE){
  result=lambda*sum(abs(vector)*abs(weights))
  return(result)
}





##loss function normal- Q
Q_normal<-function(X,y, beta, gamma_vec, delta_vec, lambda_beta, lambda_gamma, lambda_delta, w_beta, w_gamma, w_delta,l1,l2,l3,
                   already_multiplied=TRUE)
{ if (already_multiplied ==TRUE)
{gamma_vec<-array(1, dim=length(gamma_vec))
delta_vec<-array(1, dim=length(delta_vec))}
 error<-sum((y-g_normal(X=X, beta=beta, gamma_vec = gamma_vec, delta_vec = delta_vec, l1=l1, l2=l2, l3=l3, already_multiplied = already_multiplied))**2)
 penalty_beta<-get_penalty(vector=beta[unlist(get_ranges(l1,l2,l3)[1])], weights=w_beta, lambda = lambda_beta  )
 penalty_gamma<-get_penalty(vector=beta[unlist(get_ranges(l1,l2,l3)[2])]*gamma_vec, weights=w_gamma, lambda = lambda_gamma  )
 penalty_delta<-get_penalty(vector=beta[unlist(get_ranges(l1,l2,l3)[3])]*delta_vec, weights=w_delta, lambda = lambda_delta  )
 loss<- error+penalty_beta+penalty_gamma+penalty_delta
 cat("err,", error, '  ',penalty_beta,' ',penalty_gamma,' ',penalty_delta )
 return(loss)
}





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
gamma_vec_2way<-array(c(1,2,3,4,1),dim=l1*(l2+l3)+l2*l3)
delta_vec_3way<-array(c(2,3),dim=l1*l2*l3)

X=matrix(1, nrow=4, ncol=4)
X[2,4]=0
beta=array(c(1,1,1,2))
beta_vec_2way<-get_beta_vec_2way(beta=beta,l1=l1,l2=l2,l3=l3, gamma=gamma_vec_2way, only_beta = FALSE)
beta_vec_3way<-get_beta_vec_3way(beta_vec_2way,l1,l2,l3, delta=delta_vec_3way, only_beta = FALSE)
xx.all<-get_xx.all(X,l1,l2,l3)
xxx.all<-get_xxx.all(X,l1,l2,l3)
lambda_beta<-2
lambda_gamma<-3
lambda_delta<-4
w_beta<-array(2, dim=length(beta))
w_gamma<-array(3, dim=length(gamma_vec_2way))
w_gamma[1]=0
w_delta<-array(3, dim=length(delta_vec_3way))
w_delta[2]=0
y=matrix(array(c(-1,-1,0,1)), ncol = 1)


beta.all<-array(c(beta, beta_vec_2way, beta_vec_3way))

g_normal(xxx.all,beta.all, gamma_vec_2way, delta_vec_3way, l1,l2,l3, already_multiplied = TRUE)
get_penalty(beta, weights=c(1,2,3,4), lambda=0.5)

Q_normal(X=xxx.all,y=y, beta=beta.all, gamma_vec = gamma_vec_2way, delta_vec = delta_vec_3way, lambda_beta = lambda_beta, lambda_gamma = lambda_gamma,
         lambda_delta = lambda_delta, w_beta = w_beta, w_gamma = w_gamma, w_delta = w_delta,l1=l1,l2=l2,l3=l3, already_multiplied = TRUE) ##inlocuieste si pleaca cu beta simplu

## SHIM CLASS for 3 way
SHIM_3way<-function(X,y, beta_init, gamma_init, delta_init,l1=36,l2=3,l3=4, scale=FALSE)
  
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
  
  
  
  fit<-function(X, y, lambda_beta, lambda_gamma, lambda_delta, w_beta=NULL, w_gamma=NULL, w_delta=NULL,  tol=1e-6, max_iter=50, compute_Q=Q_normal)
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
      Q_new<-compute_Q(X=X,y=y, beta= beta_hat, gamma_vec=gamma_hat, delta_vec=detlta_hat, 
                       lambda_beta=lambda_beta, lambda_gamma=lambda_gamma, lambda_delta=lambda_delta,
                       w_beta =w_beta, w_gamma=w_gamma, w_delta=w_delta)
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
  
  
  

