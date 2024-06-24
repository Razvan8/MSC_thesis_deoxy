library(glmnet)
source("Create_synthetic_datasets.R")
library(lars)
library(polynom)
source("helper_functions.R")



## FUNCTIONS USED FOR THE CLASS

r2 <- function(actual, predicted) {
  # Calculate the mean of the actual values
  mean_actual <- mean(actual)
  
  # Calculate the total sum of squares
  total_sum_squares <- sum((actual - mean_actual)^2)
  
  # Calculate the residual sum of squares
  residual_sum_squares <- sum((actual - predicted)^2)
  
  # Calculate R-squared
  r_squared <- 1 - (residual_sum_squares / total_sum_squares)
  
  return(r_squared)
}

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


Soft_thresholding <- function(c, lambda) {
  assert <- function(condition, message) {
    if (!condition) stop(message)
  }
  
  assert(lambda >= 0, "lambda cannot be negative.")
  
  # Apply soft thresholding component-wise
  result <- sign(c) * pmax(abs(c) - c(lambda), 0)
  
  return(result)
}


lasso_1d_closed_form<-function(X, y, lambda, w=1, scaled=FALSE)
{ xty<-sum(X*y)
  xtx<-sum(X*X)
  if (xtx ==0)
  {return(0)}
  c<-lambda*w
  if (scaled ==TRUE)
  {c<-c*length(X)}
  #cat(" xtx: ", xtx, " xty: ", xty, " c: ", c )
  #cat(" (xty-c) / xtx:", (xty-c)/xtx, xty/xtx -c/xtx)
  
  result<-Soft_thresholding(xty, c/2)/xtx
  return(result)
}

lasso_R<-function(X,y,lambda,w)
{lmd<-lambda*w
zeros<-X*0+rnorm(length(X), 1, 0)
X<-cbind(X,zeros)
final_model <- glmnet(X, y, alpha = 1, lambda = lmd, intercept = FALSE, standardize = FALSE)
#print( coef(final_model))
coefficients <- coef(final_model)[1]
return(coefficients)}

# 2,2 4.4
# conj 3.03 5.2
#lass 2.81

##poly_minimum for 4th degree with sign



poly_min_sign<-function(coefs, positive) #coefs should be from c0 to c4
{assert(length(coefs)==5, "A 4th deg poly should have 5 coefs")
  

  c4<-coefs[5]
  c3<-coefs[4]
  c2<-coefs[3]
  c1<-coefs[2]
  c0<-coefs[1]
  coefs<-array(coefs)
  poly <- function(x) {
    x_vect<-c(1,x,x^2,x^3,x^4)
    return(sum(x_vect*coefs))
  }
  d3<-4*c4
  d2<-3*c3
  d1<-2*c2
  d0<-c1
  coefs_deriv<-c(d0,d1,d2,d3)
  #print(coefs_deriv)
  roots_deriv <- polyroot(coefs_deriv)
  real_roots<-roots_deriv[abs(Im(roots_deriv))<=1e-14]
  real_roots<-Re(real_roots)
  if (positive == TRUE)
  {if(length(real_roots[real_roots>0])==0)
  {real_roots<-list(1e-10)} #add one as backup if it does not find a root
    real_roots<-real_roots[real_roots>0]}
  if(positive == FALSE)
  {if(length(real_roots[real_roots<0])==0)
  {real_roots<-list(-1e-10)} # add one as backup if it does not find a root
    real_roots<-real_roots[real_roots<0]}
  
  
  # Evaluate the polynomial at the real roots
  values <- sapply(real_roots, function(x) poly(Re(x)))
  #print(values)
  # Find the minimum value and the corresponding root
  min_value <- min(values)
  min_root <- real_roots[which.min(values)]
  
  # Print the results
  #print("root for min:")
  #print(min_root)  # The x-value at which the minimum occurs
  #print("Min value of function")
  #print(min_value)  # The minimum value of the polynomial
  
  return(c(min_root, min_value))
  
  
  }

poly_min_sign(c(4,4,-3,-2,1), positive = FALSE)



poly_lasso_min<-function(coefs, lambda){ #finds the min for a< a=0, a>0 gets the min out of all
  ##long term take care if 2 with same min smart decision
  
  coefs_without_lambda<-coefs
  
  #case 1 sign positive
  coefs_poz<-coefs_without_lambda
  coefs_poz[2]<-coefs_poz[2]+lambda
  poz<-poly_min_sign(coefs=coefs_poz, positive = TRUE)
  x_min_poz<-poz[1]
  val_min_poz<-poz[2]
  
  #print("starts neg")

  #case 2 sign negative
  coefs_neg<-coefs_without_lambda
  coefs_neg[2]<-coefs_neg[2]-lambda
  neg<-poly_min_sign(coefs=coefs_neg, positive = FALSE)
  x_min_neg<-neg[1]
  val_min_neg<-neg[2]
  
  #print("start 0")
  
  #case 3 sign is 0
  x_zeros<-0
  val_min_zeros<-coefs[1]
  
  all_mins<-c(val_min_neg, val_min_poz, val_min_zeros)
  all_x_min<-c(x_min_neg, x_min_poz, x_zeros)#
  idx_min<-which.min(all_mins)
  x_min<-all_x_min[idx_min]
  #print("x mins neg poz 0:")
  #print(all_x_min)
  #print(" min neg poz 0: ")
  #print(all_mins)
  #print("RESULTS x_min, f min: ")
  #print(x_min)
  #print(all_mins[idx_min])
  return(unlist(x_min))
  
}


poly_lasso_min(c(1/16,1/2,3/2,2,1), lambda=1e-0)


get_coef_from_xyz<-function(x,y,z)
{c4<-sum(z*z)
 c3<-sum(2*sum(x*z))
 c2<-sum(x*x)-2*sum(y*z)
 c1<--2*sum(x*y)
 c0<-sum(y*y)
 return(c(c0,c1,c2,c3,c4))}


#get_coef_from_xyz(x=array(0, dim=10), y=array(1, dim=10), z=c(1,2,1,2,1,2,1,2,1,2))

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


get_positions_2way<-function(ls_positions, l1, l2, l3){
  all_positions<-c()
  for (tuple in ls_positions)
  {pos<-matrix_position_to_vector_index_2way(position_tuple = tuple, l1=l1, l2=l2, l3=l3)
  all_positions<-c(all_positions, pos)}
  return(all_positions)}
  
#ls_pos<-list(c(1,3), c(1,4), c(2,4), c(2,3), c(3,6))  
#get_positions_2way(ls_pos, 2,2,2)




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

get_positions_3way<-function(ls_positions, l1, l2, l3)
{all_positions<-c()
for (tuple in ls_positions)
{pos<-table_position_to_vector_index3(position_tuple = tuple, l1=l1, l2=l2, l3=l3)
#print(pos)
all_positions<-c(all_positions, pos)}
return(all_positions)}

#ls_pos<-list(c(1,3,6), c(1,3,5), c(1,4,7), c(2,4,7))
#get_positions_3way(ls_pos, l1=2, l2=2, l3=3)

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



get_beta_vec_2way<-function(beta,l1,l2,l3, gamma, only_beta = FALSE) #beta is beta_main
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


get_beta_vec_3way<-function(beta_2way,l1,l2,l3, delta, only_beta=FALSE) ## # beta_2way should be final beta_2way; #only_beta means product of beta_2ways (gamma included) without delta
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

two_ways_contribution<-function(X, gamma_vec, beta_vec_2way,l1,l2,l3, already_multiplied=FALSE) ##assumes gamma_vec is in same order as X, beta vec is only beta
{if (already_multiplied==TRUE)
{gamma_vec<-array(1, dim=length(gamma_vec))}
  range_2ways<-unlist(get_ranges(l1,l2,l3)[2])
#print(dim())
 two_ways_contrib<- X[,range_2ways]%*%(beta_vec_2way*gamma_vec) ##last multiplication should be elementwise
 return(two_ways_contrib)}

three_ways_contribution<-function(X, delta_vec, beta_vec_3way, l1,l2,l3, already_multiplied=FALSE) ##assumes gamma_vec is in same order as X and beta_vec_3way only prod of beta2way
{if (already_multiplied==TRUE)
{delta_vec<-array(1, dim=length(delta_vec))}
  range_3ways<-unlist(get_ranges(l1,l2,l3)[3])
three_ways_contrib<-X[,range_3ways]%*%(beta_vec_3way*delta_vec) ##last multiplication should be elementwise
return(three_ways_contrib)}



### g function

g_normal<-function(X, beta, gamma_vec, delta_vec,l1,l2,l3, already_multiplied=TRUE) #bet_2way is without gamma, beta_3way withour delta only
  #already multiplied=True means beta already has gamma delta inside
{beta_main<-beta[unlist(get_ranges(l1,l2,l3)[1])]
 beta_2way<-beta[unlist(get_ranges(l1,l2,l3)[2])]
 beta_3way<-beta[unlist(get_ranges(l1,l2,l3)[3])]
 if (already_multiplied==TRUE)
 {gamma_vec<-array(1, dim=length(gamma_vec))
  delta_vec<-array(1, dim=length(delta_vec))}
 result<-mains_contribution(X=X, beta_main = beta_main, l1=l1, l2=l2, l3=l3)+ 
         two_ways_contribution(X=X, gamma_vec=gamma_vec, beta_vec_2way=beta_2way,l1=l1, l2=l2, l3=l3, already_multiplied = FALSE)+
         three_ways_contribution(X=X, delta_vec = delta_vec, beta_vec_3way = beta_3way,l1=l1, l2=l2, l3=l3, already_multiplied = FALSE)
 #cat("g:", result[1:10])
 return(result)
}



##penalty for 1 vector
get_penalty<-function(vector, weights, lambda, already_weighted=TRUE){
  result=lambda*sum(abs(vector)*abs(weights))
  return(result)
}





##loss function normal- Q
Q_normal<-function(X,y, beta, gamma_vec, delta_vec, lambda_beta, lambda_gamma, lambda_delta, w_beta, w_gamma, w_delta,l1,l2,l3,
                   already_multiplied=TRUE, scaled=FALSE)
{ if (length(beta)== l1 +l2+l3)
{#print("Beta was given only main and computed for the rest")
  already_multiplied = TRUE
  beta_2way<-get_beta_vec_2way(beta = beta, l1=l1, l2=l2, l3=l3, gamma= gamma_vec, only_beta = FALSE )
  beta_3way<-get_beta_vec_3way(beta = beta_2way, l1=l1, l2=l2, l3=l3, delta = delta_vec, only_beta = FALSE)
  beta<-c(beta, beta_2way,beta_3way)}
  #should be before making it 1
  penalty_beta<-get_penalty(vector=beta[unlist(get_ranges(l1,l2,l3)[1])], weights=w_beta, lambda = lambda_beta  )
  penalty_gamma<-get_penalty(vector=gamma_vec, weights=w_gamma, lambda = lambda_gamma  )
  penalty_delta<-get_penalty(vector=delta_vec, weights=w_delta, lambda = lambda_delta  )
  
  if (already_multiplied ==TRUE)
{gamma_vec<-array(1, dim=length(gamma_vec))
delta_vec<-array(1, dim=length(delta_vec))}
 error<-sum((y-g_normal(X=X, beta=beta, gamma_vec = gamma_vec, delta_vec = delta_vec, l1=l1, l2=l2, l3=l3, already_multiplied = already_multiplied))**2)
  if(scaled==TRUE)
   {error<-error/(2*dim(X)[1])}
 loss<- error+penalty_beta+penalty_gamma+penalty_delta
 #cat("err,", error, '  ',penalty_beta,' ',penalty_gamma,' ',penalty_delta )
 return(loss)
}





###RELATIVE DIFFERENCE

compute_rel_dif<-function(Q_old, Q_new)
{rel_dif<- abs(Q_old-Q_new)/abs(Q_old)
return(rel_dif)}



##### UPDATE DELTA FUNCTION #####

update_delta<-function(X, y,beta_hat, gamma_hat, delta_hat, lambda_delta, l1, l2, l3) 
  ##beta_hat is only for mains
{beta_2way <- get_beta_vec_2way(beta = beta_hat, l1=l1, l2=l2, l3=l3, gamma=gamma_hat, only_beta = FALSE) ###This is with delta
 beta_3way <- get_beta_vec_3way(beta_2way = beta_2way, l1=l1, l2=l2, l3=l3, delta=delta_hat, only_beta = TRUE) #This is with gamma WITHOUTdelta
 y_tilde <- y - mains_contribution(X=X, beta_main = beta_hat, l1=l1,l2=l2, l3=l3) 
           -two_ways_contribution(X=X, gamma_vec = gamma_hat, beta_vec_2way = beta_2way, l1=l1, l2=l2, l3=l3, already_multiplied = TRUE )
 X_3way<-X[,unlist(get_ranges(l1,l2,l3)[3])]

 X_tilde<-matrix(rep(beta_3way, each = nrow(X_3way)), nrow = nrow(X_3way))*X_3way
 lambda_delta<-lambda_delta/(2*nrow(X)) ##scale lambda because in lasso we have 1/(2n)* loss
 lasso_model <- glmnet(X_tilde, y_tilde, alpha = 1, lambda = lambda_delta, intercept = FALSE, standardize = FALSE)
 #print(lambda_delta)

 lasso_coef <- coef(lasso_model)
 delta_hat<- lasso_coef[-1]
 return(delta_hat)
 
}



###test update delta
data<- create_basic_dataset()
X<- data$X
y<- data$y$true

beta_true<- data$beta[-1,]
l1=8
l2=8
l3=4

for (lambda_delta in c( c(0), c(100) ) ) {
  for (noise in c(0,0.01,1)){

beta_main<-beta_true[1:(l1+l2+l3)]
beta_2way<-beta_true[unlist(get_ranges(l1,l2,l3)[2])]
beta_3way<-beta_true[unlist(get_ranges(l1,l2,l3)[3])]
beta_2way_without_gamma<-get_beta_vec_2way(beta = beta_main, l1=l1, l2=l2, l3=l3, gamma=NULL, only_beta = TRUE )
gamma_hat<-beta_2way/beta_2way_without_gamma
gamma_hat[is.nan(gamma_hat)]<-0
#lambda_delta<-100

beta_3way_without_gamma<-get_beta_vec_3way(beta_2way = beta_2way, l1=l1, l2=l2, l3=l3, delta = NULL, only_beta = TRUE)

delta_true<-beta_3way/beta_3way_without_gamma
delta_true[is.nan(delta_true)]<-0

delta_hat<-delta_true+rnorm(length(delta_true),0,noise)
delta_pred <- update_delta(X=X, y=y,beta_hat=beta_main, gamma_hat=gamma_hat, delta_hat=delta_hat, lambda_delta=lambda_delta, l1=l1, l2=l2, l3=l3) 
 


beta_3way_with_gamma<-get_beta_vec_3way(beta_2way = beta_2way, l1=l1, l2=l2, l3=l3, delta = delta_hat, only_beta = FALSE)




q_pred <- Q_normal(X=X,y=y, beta=beta_main, gamma_vec=gamma_hat, delta_vec=delta_pred, 
         lambda_beta=1, lambda_gamma=1, lambda_delta=lambda_delta, 
          w_beta=1, w_gamma=1, w_delta=1,l1=l1,l2=l2,l3=l3, already_multiplied=TRUE)

q_before<-Q_normal(X=X,y=y, beta=beta_main, gamma_vec=gamma_hat, delta_vec=delta_hat, 
         lambda_beta=1, lambda_gamma=1, lambda_delta=lambda_delta, 
         w_beta=1, w_gamma=1, w_delta=1,l1=l1,l2=l2,l3=l3,already_multiplied=TRUE)

print("info: ")
cat(" noise:",noise,"lmd:", lambda_delta)
print(q_pred-q_before)

if(q_pred-q_before>0)
{print("numerical instability !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")}

#print(delta_pred)
#print("pred was now true")
#print(delta_true)

}}

print(delta_pred)



##### UPDATE GAMMA FUNCTION #####


update_gamma<-function(X, y,beta_hat, gamma_hat, delta_hat, lambda_gamma, l1, l2, l3, w=1) 
{range1<- c(1:l1)
 range2<-c((l1+1):(l1+l2))
 range3<-c( (l1+l2+1) : (l1+l2+l3) )
 X_2way<-X[,c( (l1+l2+l3+1): (l1+l2+l3+l1*l2+l2*l3+l1*l3) )]
 X_3way<-X[,c( (l1+l2+l3+l1*l2+l2*l3+l1*l3+1):(l1+l2+l3+l1*l2+l2*l3+l1*l3+l1*l2*l3) )]
 beta_2way <- get_beta_vec_2way(beta = beta_hat, l1=l1, l2=l2, l3=l3, gamma=gamma_hat, only_beta = FALSE) ###This is with delta
 beta_3way <- get_beta_vec_3way(beta_2way = beta_2way, l1=l1, l2=l2, l3=l3, delta=delta_hat, only_beta = FALSE) #This is with gamma WITH delta
 
 if (w==1)
 {w=array(1, dim=length(gamma_hat))}
 
for(i in range1){
  for (j in range2){
    
    beta_2way <- get_beta_vec_2way(beta = beta_hat, l1=l1, l2=l2, l3=l3, gamma=gamma_hat, only_beta = FALSE) ###This is with delta
    beta_3way <- get_beta_vec_3way(beta_2way = beta_2way, l1=l1, l2=l2, l3=l3, delta=delta_hat, only_beta = FALSE) #This is with gamma WITH delta
    
    
    discard_2way<-matrix_position_to_vector_index_2way(c(i,j), l1=l1, l2=l2, l3=l3)
    ls_tuples_discard_3way<-list()
    
    for (elem3 in range3) {
      tuple <- c(i, j, elem3)
      ls_tuples_discard_3way <- append(ls_tuples_discard_3way, list(tuple))
    }
    #print("ls tuples: ")
    #print(ls_tuples_discard_3way)
    
    discard_3way<- get_positions_3way(ls_tuples_discard_3way, l1=l1, l2=l2, l3=l3) #positions 3 way in vector form
    
    #print("ls positions 3 way: ")
    #print(discard_3way)
    X_2way_kept<-X_2way[,-discard_2way]
    X_3way_kept<-X_3way[, -discard_3way]
    beta_2way_kept <- beta_2way[-discard_2way]
    beta_3way_kept <- beta_3way[-discard_3way]
    gamma_hat_kept <- gamma_hat[-discard_2way]
    delta_hat_kept <- delta_hat[-discard_3way]
    #print(X_2way_kept)
    #print(beta_2way_kept)
    #print(X_3way_kept)
    #print(beta_3way_kept)
    
    y_tilde<-y - mains_contribution(X=X, beta_main=beta_hat, l1=l1, l2=l2, l3=l3) - X_2way_kept%*%beta_2way_kept -X_3way_kept%*%beta_3way_kept
    #print(y_tilde)
    three_ways=0
    #print("ok")
    for (k in range3) #compute 3 ways contrib
    {three_ways<-three_ways+ X_3way[,table_position_to_vector_index3(c(i,j,k),l1=l1, l2=l2, l3=l3)]*((beta_hat[i]*beta_hat[j]*beta_hat[k])^2)*
                             gamma_hat[matrix_position_to_vector_index_2way(c(i,k), l1=l1, l2=l2 ,l3=l3)] *  
                             gamma_hat[matrix_position_to_vector_index_2way(c(j,k), l1=l1, l2=l2 ,l3=l3)] *
                             delta_hat[table_position_to_vector_index3(c(i,j,k), l1=l1, l2=l2, l3=l3)]
    
    }
 
    X_tilde<-X_2way[,discard_2way]*beta_hat[i]*beta_hat[j]+  three_ways
    
    Q_old <- Q_normal(X=X,y=y, beta=beta_hat, gamma_vec=gamma_hat, delta_vec=delta_hat, 
                  lambda_beta=1, lambda_gamma=lambda_gamma, lambda_delta=1, 
                  w_beta=1, w_gamma=1, w_delta=1,l1=l1,l2=l2,l3=l3, already_multiplied=TRUE)
    
    
    gamma_hat[matrix_position_to_vector_index_2way(c(i,j), l1=l1, l2=l2, l3=l3)]<-lasso_1d_closed_form(X=X_tilde, y= y_tilde,
                                       lambda=lambda_gamma, w=w[matrix_position_to_vector_index_2way(c(i,j), l1=l1, l2=l2, l3=l3) ] )
    #gamma_hat[matrix_position_to_vector_index_2way(c(i,j), l1=l1, l2=l2, l3=l3)]<-lasso_R(X=X_tilde, y=y_tilde, lambda=lambda_gamma, w=1)
    

    
   
    
    
    beta_2way <- get_beta_vec_2way(beta = beta_hat, l1=l1, l2=l2, l3=l3, gamma=gamma_hat, only_beta = FALSE) ###This is with delta
    beta_3way <- get_beta_vec_3way(beta_2way = beta_2way, l1=l1, l2=l2, l3=l3, delta=delta_hat, only_beta = FALSE) #This is with gamma WITH delta
    
    Q_new <- Q_normal(X=X,y=y, beta=beta_hat, gamma_vec=gamma_hat, delta_vec=delta_hat, 
             lambda_beta=1, lambda_gamma=lambda_gamma, lambda_delta=1, 
             w_beta=1, w_gamma=1, w_delta=1,l1=l1,l2=l2,l3=l3,already_multiplied=TRUE)
    
    #if (Q_new-Q_old >=0)
    #cat(" new-old: ",Q_new-Q_old, " Q: ",Q_new)
    if (Q_new-Q_old >= Q_old/100)
    {print("There might be numerical instability in gamma.")}
    
    
  }
}
  

 for(i in range1){
   for (k in range3){
     
     beta_2way <- get_beta_vec_2way(beta = beta_hat, l1=l1, l2=l2, l3=l3, gamma=gamma_hat, only_beta = FALSE) ###This is with delta
     beta_3way <- get_beta_vec_3way(beta_2way = beta_2way, l1=l1, l2=l2, l3=l3, delta=delta_hat, only_beta = FALSE) #This is with gamma WITH delta
     
     discard_2way<-matrix_position_to_vector_index_2way(c(i,k), l1=l1, l2=l2, l3=l3)
     ls_tuples_discard_3way<-list()
     
     for (elem3 in range2) {
       tuple <- c(i, elem3, k)
       ls_tuples_discard_3way <- append(ls_tuples_discard_3way, list(tuple))
     }
     #print("ls tuples: ")
     #print(ls_tuples_discard_3way)
     
     discard_3way<- get_positions_3way(ls_tuples_discard_3way, l1=l1, l2=l2, l3=l3) #positions 3 way in vector form
     
     #print("ls positions 3 way: ")
     #print(discard_3way)
     X_2way_kept<-X_2way[,-discard_2way]
     X_3way_kept<-X_3way[, -discard_3way]
     beta_2way_kept <- beta_2way[-discard_2way]
     beta_3way_kept <- beta_3way[-discard_3way]
     gamma_hat_kept <- gamma_hat[-discard_2way]
     delta_hat_kept <- delta_hat[-discard_3way]
     
     y_tilde<-y - mains_contribution(X=X, beta_main=beta_hat, l1=l1, l2=l2, l3=l3) - X_2way_kept%*%beta_2way_kept -X_3way_kept%*%beta_3way_kept
     three_ways=0
     for (j in range2) #compute 3 ways contrib
     {three_ways<-three_ways+ X_3way[,table_position_to_vector_index3(c(i,j,k),l1=l1, l2=l2, l3=l3)]*((beta_hat[i]*beta_hat[j]*beta_hat[k])^2)*
       gamma_hat[matrix_position_to_vector_index_2way(c(i,j), l1=l1, l2=l2 ,l3=l3)] *  
       gamma_hat[matrix_position_to_vector_index_2way(c(j,k), l1=l1, l2=l2 ,l3=l3)] *
       delta_hat[table_position_to_vector_index3(c(i,j,k), l1=l1, l2=l2, l3=l3)]
     
     }
     X_tilde<-X_2way[,discard_2way]*beta_hat[i]*beta_hat[k]+ three_ways
     
     Q_old <- Q_normal(X=X,y=y, beta=beta_hat, gamma_vec=gamma_hat, delta_vec=delta_hat, 
                       lambda_beta=1, lambda_gamma=lambda_gamma, lambda_delta=1, 
                       w_beta=1, w_gamma=1, w_delta=1,l1=l1,l2=l2,l3=l3, already_multiplied=TRUE)
     
     
     lasso_1d_closed_form(X=X_tilde, y= y_tilde, lambda=lambda_gamma, w=w[matrix_position_to_vector_index_2way(c(i,k), l1=l1, l2=l2, l3=l3)] )
     #print(lambda_gamma)
     gamma_hat[matrix_position_to_vector_index_2way(c(i,k), l1=l1, l2=l2, l3=l3)]<-lasso_1d_closed_form(X=X_tilde, y= y_tilde, lambda=lambda_gamma, w=w[matrix_position_to_vector_index_2way(c(i,k), l1=l1, l2=l2, l3=l3) ] )
     #gamma_hat[matrix_position_to_vector_index_2way(c(i,k), l1=l1, l2=l2, l3=l3)]<-lasso_R(X=X_tilde, y=y_tilde, lambda=lambda_gamma, w=1)

     beta_2way <- get_beta_vec_2way(beta = beta_hat, l1=l1, l2=l2, l3=l3, gamma=gamma_hat, only_beta = FALSE) ###This is with delta
     beta_3way <- get_beta_vec_3way(beta_2way = beta_2way, l1=l1, l2=l2, l3=l3, delta=delta_hat, only_beta = FALSE) #This is with gamma WITH delta
     
     Q_new <- Q_normal(X=X,y=y, beta=beta_hat, gamma_vec=gamma_hat, delta_vec=delta_hat, 
                       lambda_beta=1, lambda_gamma=lambda_gamma, lambda_delta=1, 
                       w_beta=1, w_gamma=1, w_delta=1,l1=l1,l2=l2,l3=l3, already_multiplied=TRUE)
     
     #if (Q_new-Q_old >=0)
     #cat(" new-old: ",Q_new-Q_old, " Q: ",Q_new)
     if (Q_new-Q_old >= Q_old/100)
     {print("There might be numerical instability in gamma.")}
     
     
   }
 }
 


 
 for(j in range2){
   for (k in range3){
     discard_2way<-matrix_position_to_vector_index_2way(c(j,k), l1=l1, l2=l2, l3=l3)
     ls_tuples_discard_3way<-list()
     
     beta_2way <- get_beta_vec_2way(beta = beta_hat, l1=l1, l2=l2, l3=l3, gamma=gamma_hat, only_beta = FALSE) ###This is with delta
     beta_3way <- get_beta_vec_3way(beta_2way = beta_2way, l1=l1, l2=l2, l3=l3, delta=delta_hat, only_beta = FALSE) #This is with gamma WITH delta
     
     
     for (elem3 in range1) {
       tuple <- c(elem3, j, k)
       ls_tuples_discard_3way <- append(ls_tuples_discard_3way, list(tuple))
     }
     #print("ls tuples: ")
     #print(ls_tuples_discard_3way)
     
     discard_3way<- get_positions_3way(ls_tuples_discard_3way, l1=l1, l2=l2, l3=l3) #positions 3 way in vector form
     
     #print("ls positions 3 way: ")
     #print(discard_3way)
     X_2way_kept<-X_2way[,-discard_2way]
     X_3way_kept<-X_3way[, -discard_3way]
     beta_2way_kept <- beta_2way[-discard_2way]
     beta_3way_kept <- beta_3way[-discard_3way]
     gamma_hat_kept <- gamma_hat[-discard_2way]
     delta_hat_kept <- delta_hat[-discard_3way]
     
     y_tilde<-y - mains_contribution(X=X, beta_main=beta_hat, l1=l1, l2=l2, l3=l3) - X_2way_kept%*%beta_2way_kept -X_3way_kept%*%beta_3way_kept
     three_ways=0
     for (i in range1) #compute 3 ways contrib
     {three_ways<-three_ways+ X_3way[,table_position_to_vector_index3(c(i,j,k),l1=l1, l2=l2, l3=l3)]*((beta_hat[i]*beta_hat[j]*beta_hat[k])^2)*
       gamma_hat[matrix_position_to_vector_index_2way(c(i,j), l1=l1, l2=l2 ,l3=l3)] *  
       gamma_hat[matrix_position_to_vector_index_2way(c(i,k), l1=l1, l2=l2 ,l3=l3)] *
       delta_hat[table_position_to_vector_index3(c(i,j,k), l1=l1, l2=l2, l3=l3)]
     
     }
     X_tilde<-X_2way[,discard_2way]*beta_hat[j]*beta_hat[k]+ three_ways
     
     Q_old <- Q_normal(X=X,y=y, beta=beta_hat, gamma_vec=gamma_hat, delta_vec=delta_hat, 
                       lambda_beta=1, lambda_gamma=lambda_gamma, lambda_delta=1, 
                       w_beta=1, w_gamma=1, w_delta=1,l1=l1,l2=l2,l3=l3, already_multiplied=TRUE)
     
     
     gamma_hat[matrix_position_to_vector_index_2way(c(j,k), l1=l1, l2=l2, l3=l3)]<- lasso_1d_closed_form(X=X_tilde, y= y_tilde,
                                                                                                         lambda= lambda_gamma, w=w[matrix_position_to_vector_index_2way(c(j,k), l1=l1, l2=l2, l3=l3) ] )

     beta_2way <- get_beta_vec_2way(beta = beta_hat, l1=l1, l2=l2, l3=l3, gamma=gamma_hat, only_beta = FALSE) ###This is with delta
     beta_3way <- get_beta_vec_3way(beta_2way = beta_2way, l1=l1, l2=l2, l3=l3, delta=delta_hat, only_beta = FALSE) #This is with gamma WITH delta
     
     Q_new <- Q_normal(X=X,y=y, beta=beta_hat, gamma_vec=gamma_hat, delta_vec=delta_hat, 
                       lambda_beta=1, lambda_gamma=lambda_gamma, lambda_delta=1, 
                       w_beta=1, w_gamma=1, w_delta=1,l1=l1,l2=l2,l3=l3, already_multiplied=TRUE)
     #if (Q_new-Q_old >=0)
     #cat(" new-old: ",Q_new-Q_old, " Q: ",Q_new)
     if (Q_new-Q_old >= Q_old/100)
     {print("There might be numerical instability in gamma.")}
     
     
   }
 }
 
return(gamma_hat)
  
}


  
data<- create_basic_dataset()
X<- data$X
y<- data$y$true

beta_true<- data$beta[-1,]
l1=8
l2=8
l3=4
print(beta_true)

beta_main<-beta_true[1:(l1+l2+l3)]
beta_2way<-beta_true[unlist(get_ranges(l1,l2,l3)[2])]
beta_3way<-beta_true[unlist(get_ranges(l1,l2,l3)[3])]
beta_2way_without_gamma<-get_beta_vec_2way(beta = beta_main, l1=l1, l2=l2, l3=l3, gamma=NULL, only_beta = TRUE )
beta_3way_without_gamma<-get_beta_vec_3way(beta_2way = beta_2way, l1=l1, l2=l2, l3=l3, delta = NULL, only_beta = TRUE)

gamma_hat<-beta_2way/beta_2way_without_gamma
gamma_hat[is.nan(gamma_hat)]<-0
lambda_gamma<-1


gamma_true<-beta_2way/beta_2way_without_gamma
gamma_true[is.nan(gamma_true)]<-0

delta_true<-beta_3way/beta_3way_without_gamma
delta_true[is.nan(delta_true)]<-0

delta_hat<-delta_true

lambda_gamma<-2e3

gamma_hat<-gamma_true+rnorm(length(gamma_hat), 0,0.1)
gamma_hat[1]<-30
gamma_hat
gamma_pred <- update_gamma(X=X, y=y,beta_hat=beta_main, gamma_hat=gamma_hat, delta_hat=delta_hat, lambda_gamma=lambda_gamma,  l1=l1, l2=l2, l3=l3, w=1) 
gamma_pred

sum(gamma_pred==0)  
sum(abs(gamma_pred))


sum(abs(gamma_true))
sum(gamma_true==0)




##### UPDATE BETA FUNCTION #####


update_beta <- function(X, y, beta_hat, gamma_hat, delta_hat, lambda_beta, l1, l2, l3, w=1)
{
  range1 <- c(1:l1)
  range2 <- c((l1 + 1):(l1 + l2))
  range3 <- c((l1 + l2 + 1):(l1 + l2 + l3))
  X_main<-  X[, c( 1:(l1 + l2 + l3 ) )]
  X_2way <- X[, c((l1 + l2 + l3 + 1):(l1 + l2 + l3 + l1 * l2 + l2 * l3 + l1 * l3))]
  X_3way <-X[, c((l1 + l2 + l3 + l1 * l2 + l2 * l3 + l1 * l3 + 1):(l1 + l2 + l3 + l1 *l2 + l2 * l3 + l1 * l3 + l1 * l2 * l3))]
  beta_2way <- get_beta_vec_2way(beta = beta_hat, l1 = l1, l2 = l2, l3 = l3, gamma = gamma_hat, only_beta = FALSE) ###This is with delta
  beta_3way <- get_beta_vec_3way(beta_2way = beta_2way, l1 = l1, l2 = l2, l3 = l3, delta = delta_hat, only_beta = FALSE) #This is with gamma WITH delta
  
  if (w == 1)
  {
    w = array(1, dim = length(gamma_hat))
  }
  
  ##CASE 1 Iterate for i
  for (i in range1) {
    discard_main <- c(i)
    #print(i)
    
    #discard 2way
    ls_pos_2way <- list()
    for (jk in c(range2, range3))
      #positions ij ik in this order
    {
      ls_pos_2way <- append(ls_pos_2way, list(c(i, jk)) )}
    #print("before get_poss")
    #print(ls_pos_2way)
    discard_2way <- get_positions_2way(ls_positions = ls_pos_2way, l1 = l1, l2 = l2,l3 = l3) #discard 2 way in vector indexing
    #print("after get_poss")
    
    
    #discard 3way
    ls_pos_3way<-list()
    for (j in range2) {
      for (k in range3) {
        ls_pos_3way <- append(ls_pos_3way, list(c(i, j, k)) )
        
      }
    }
    discard_3way <- get_positions_3way(ls_positions = ls_pos_3way,l1 = l1, l2 = l2, l3 = l3)
    
    #print("X3way")
    #print(X_3way)
    #print("dicard_3way")
    #print(discard_3way)
    
    X_main_kept<-X_main[,-discard_main]
    X_2way_kept <- X_2way[, -discard_2way]
    X_3way_kept <- X_3way[,-discard_3way]
    
    beta_main_kept<-beta_hat[-discard_main]
    beta_2way_kept <- beta_2way[-discard_2way]
    beta_3way_kept <- beta_3way[-discard_3way]
    gamma_hat_kept <- gamma_hat[-discard_2way]
    delta_hat_kept <- delta_hat[-discard_3way]
    
    
    y_tilde<-y - X_main_kept%*%beta_main_kept - X_2way_kept%*%beta_2way_kept -X_3way_kept%*%beta_3way_kept
    #print("ytilde")
    #print(y_tilde)
    
    #Create X_tilde from X1_tilde X2_tilde
    two_ways<-0
    for (jk in c(range2, range3)) #compute 3 ways contrib
    {two_ways<-two_ways + X_2way[,matrix_position_to_vector_index_2way(c(i,jk),l1=l1, l2=l2, l3=l3)]*(beta_hat[jk])*
      gamma_hat[matrix_position_to_vector_index_2way(c(i,jk), l1=l1, l2=l2 ,l3=l3)] }
    
    three_ways<-0
    for (j in range2) #compute 3 ways contrib
    {for(k in range3)
    { #print("X_3way[cijk]")
      #print(X_3way[,table_position_to_vector_index3(c(i,j,k),l1=l1, l2=l2, l3=l3)])
      #print('gamma')
      #print(gamma_hat[matrix_position_to_vector_index_2way(c(i,j), l1=l1, l2=l2 ,l3=l3)])
      #print("delta")
      #print(delta_hat[table_position_to_vector_index3(c(i,j,k), l1=l1, l2=l2, l3=l3)])
      
      
      
      
      
      three_ways<-three_ways+ X_3way[,table_position_to_vector_index3(c(i,j,k),l1=l1, l2=l2, l3=l3)]*((beta_hat[j]*beta_hat[k])^2)*
      gamma_hat[matrix_position_to_vector_index_2way(c(i,j), l1=l1, l2=l2 ,l3=l3)] *  
      gamma_hat[matrix_position_to_vector_index_2way(c(i,k), l1=l1, l2=l2 ,l3=l3)] *
      gamma_hat[matrix_position_to_vector_index_2way(c(j,k), l1=l1, l2=l2 ,l3=l3)] *
      delta_hat[table_position_to_vector_index3(c(i,j,k), l1=l1, l2=l2, l3=l3)]
    
    }}
    #cat("X[, discard main]: ", X[,discard_main], "dim2ways", two_ways)
    #print(class(X_main[,discard_main]))
    #print(class(two_ways))
    
  
    
    X1_tilde<- array( X_main[,discard_main]) + array(two_ways)
    X2_tilde<-three_ways
    #print("3ways")
    #print(three_ways)
    
    #cat("X1dim: ", dim(X1_tilde), " X2dim: ", dim(X2_tilde))
    Q_old <- Q_normal(X=X,y=y, beta=beta_hat, gamma_vec=gamma_hat, delta_vec=delta_hat, 
                      lambda_beta=lambda_beta, lambda_gamma=0, lambda_delta=0, 
                      w_beta=1, w_gamma=1, w_delta=1,l1=l1,l2=l2,l3=l3, already_multiplied=TRUE)
    

    
    coefs<-get_coef_from_xyz(x=array(X1_tilde), y= array(y_tilde), z=array(X2_tilde) ) #coefs as c0 c1 c2...c4
    beta_hat[i]<-poly_lasso_min(coefs = coefs, lambda = lambda_beta) #beta updated

    beta_2way <- get_beta_vec_2way(beta = beta_hat, l1 = l1, l2 = l2, l3 = l3, gamma = gamma_hat, only_beta = FALSE) ###This is with delta
    beta_3way <- get_beta_vec_3way(beta_2way = beta_2way, l1 = l1, l2 = l2, l3 = l3, delta = delta_hat, only_beta = FALSE) #This is with gamma WITH delta
    
    
    
    Q_new <- Q_normal(X=X,y=y, beta=beta_hat, gamma_vec=gamma_hat, delta_vec=delta_hat, 
                      lambda_beta=lambda_beta, lambda_gamma=0, lambda_delta=0, 
                      w_beta=1, w_gamma=1, w_delta=1,l1=l1,l2=l2,l3=l3, already_multiplied=TRUE)
    if (Q_new-Q_old >= Q_old/100)
    {print("There might be numerical instability in beta")}

    #cat(" new-old: ",Q_new-Q_old, " Q: ",Q_new)
    
  }
  

  
  #print("start j")
  
  ##CASE 2 Iterate for j
  for (j in range2) {
    discard_main <- c(j)
    
    #discard 2way
    ls_pos_2way <- list()
    for (i in range1) #ij
    {ls_pos_2way <- append(ls_pos_2way, list(c(i,j)) )}
    for (k in range3) #ij
    {ls_pos_2way <- append(ls_pos_2way, list( c(j,k)) )}
    
    discard_2way <- get_positions_2way(ls_positions = ls_pos_2way, l1 = l1, l2 = l2,l3 = l3) #discard 2 way in vector indexing
    
    
    #discard 3way
    ls_pos_3way<-list()
    for (i in range1) {
      for (k in range3) {
        ls_pos_3way <- append(ls_pos_3way, list(c(i, j, k) ))
        
      }
    }
    #print(ls_pos_3way)
    discard_3way <- get_positions_3way(ls_positions = ls_pos_3way,l1 = l1, l2 = l2, l3 = l3)
    
    X_main_kept<-X_main[,-discard_main]
    X_2way_kept <- X_2way[, -discard_2way]
    X_3way_kept <- X_3way[,-discard_3way]
    
    beta_main_kept<-beta_hat[-discard_main]
    beta_2way_kept <- beta_2way[-discard_2way]
    beta_3way_kept <- beta_3way[-discard_3way]
    gamma_hat_kept <- gamma_hat[-discard_2way]
    delta_hat_kept <- delta_hat[-discard_3way]
    y_tilde<-y - X_main_kept%*%beta_main_kept - X_2way_kept%*%beta_2way_kept -X_3way_kept%*%beta_3way_kept
    

    
    #Create X_tilde from X1_tilde X2_tilde
    two_ways<-0
    for (i in range1) #compute 3 ways contrib
    {two_ways<-two_ways + X_2way[,matrix_position_to_vector_index_2way(c(i,j),l1=l1, l2=l2, l3=l3)]*(beta_hat[i])*
      gamma_hat[matrix_position_to_vector_index_2way(c(i,j), l1=l1, l2=l2 ,l3=l3)] }
    for (k in range3) #compute 3 ways contrib
    {two_ways<-two_ways + X_2way[,matrix_position_to_vector_index_2way(c(j,k),l1=l1, l2=l2, l3=l3)]*(beta_hat[k])*
      gamma_hat[matrix_position_to_vector_index_2way(c(j,k), l1=l1, l2=l2 ,l3=l3)] }
   
    
    three_ways<-0
    for (i in range1) #compute 3 ways contrib
    {for(k in range3)
    {three_ways<-three_ways+ X_3way[,table_position_to_vector_index3(c(i,j,k),l1=l1, l2=l2, l3=l3)]*((beta_hat[i]*beta_hat[k])^2)*
      gamma_hat[matrix_position_to_vector_index_2way(c(i,j), l1=l1, l2=l2 ,l3=l3)] *  
      gamma_hat[matrix_position_to_vector_index_2way(c(i,k), l1=l1, l2=l2 ,l3=l3)] *
      gamma_hat[matrix_position_to_vector_index_2way(c(j,k), l1=l1, l2=l2 ,l3=l3)] *
      delta_hat[table_position_to_vector_index3(c(i,j,k), l1=l1, l2=l2, l3=l3)]
    
    }}
    
    X1_tilde<- X[,discard_main] + two_ways
    X2_tilde<-three_ways
    
    
    Q_old <- Q_normal(X=X,y=y, beta=beta_hat, gamma_vec=gamma_hat, delta_vec=delta_hat, 
                      lambda_beta=lambda_beta, lambda_gamma=0, lambda_delta=0, 
                      w_beta=1, w_gamma=1, w_delta=1,l1=l1,l2=l2,l3=l3, already_multiplied=TRUE)

    
    coefs<-get_coef_from_xyz(x=array(X1_tilde), y= array(y_tilde), z=array(X2_tilde)) #coefs as c0 c1 c2...c4
    beta_hat[j]<-poly_lasso_min(coefs = coefs, lambda = lambda_beta) #beta updated

    
    #update beta_23way
    beta_2way <- get_beta_vec_2way(beta = beta_hat, l1 = l1, l2 = l2, l3 = l3, gamma = gamma_hat, only_beta = FALSE) ###This is with delta
    beta_3way <- get_beta_vec_3way(beta_2way = beta_2way, l1 = l1, l2 = l2, l3 = l3, delta = delta_hat, only_beta = FALSE) #This is with gamma WITH delta
    
    
    Q_new <- Q_normal(X=X,y=y, beta=beta_hat, gamma_vec=gamma_hat, delta_vec=delta_hat, 
                      lambda_beta=lambda_beta, lambda_gamma=0, lambda_delta=0, 
                      w_beta=1, w_gamma=1, w_delta=1,l1=l1,l2=l2,l3=l3, already_multiplied=TRUE)
    if (Q_new-Q_old >= Q_old/100)
    {print("There might be numerical instability in beta.")}
    
    #cat(" new-old: ",Q_new-Q_old, " Q: ",Q_new) 
  }
    
    
    
    ##CASE 3 Iterate for k
    for (k in range3) {
      discard_main <- c(k)
      
      #discard 2way
      ls_pos_2way <- list()
      for (ij in c(range1, range2))
        #positions ij ik in this order
      {
        ls_pos_2way <- append(ls_pos_2way, list( c(ij, k)) )}
      discard_2way <- get_positions_2way(ls_positions = ls_pos_2way, l1 = l1, l2 = l2,l3 = l3) #discard 2 way in vector indexing
      
      
      #discard 3way
      ls_pos_3way<-list()
      for (i in range1) {
        for (j in range2) {
          ls_pos_3way <- append(ls_pos_3way, list( c(i, j, k)) )
          
        }
      }
      discard_3way <- get_positions_3way(ls_positions = ls_pos_3way,l1 = l1, l2 = l2, l3 = l3)
      
      X_main_kept<-X_main[,-discard_main]
      X_2way_kept <- X_2way[, -discard_2way]
      X_3way_kept <- X_3way[,-discard_3way]
      
      beta_main_kept<-beta_hat[-discard_main]
      beta_2way_kept <- beta_2way[-discard_2way]
      beta_3way_kept <- beta_3way[-discard_3way]
      gamma_hat_kept <- gamma_hat[-discard_2way]
      delta_hat_kept <- delta_hat[-discard_3way]
      y_tilde<-y -X_main_kept%*%beta_main_kept - X_2way_kept%*%beta_2way_kept -X_3way_kept%*%beta_3way_kept
      
      #Create X_tilde from X1_tilde X2_tilde
      two_ways<-0
      for (ij in c(range1, range2)) #compute 3 ways contrib
      {two_ways<-two_ways + X_2way[,matrix_position_to_vector_index_2way(c(ij,k),l1=l1, l2=l2, l3=l3)]*(beta_hat[ij])*
        gamma_hat[matrix_position_to_vector_index_2way(c(ij,k), l1=l1, l2=l2 ,l3=l3)] }
      
      three_ways<-0
      for (i in range1) #compute 3 ways contrib
      {for(j in range2)
      {three_ways<-three_ways+ X_3way[,table_position_to_vector_index3(c(i,j,k),l1=l1, l2=l2, l3=l3)]*((beta_hat[i]*beta_hat[j])^2)*
        gamma_hat[matrix_position_to_vector_index_2way(c(i,j), l1=l1, l2=l2 ,l3=l3)] *  
        gamma_hat[matrix_position_to_vector_index_2way(c(i,k), l1=l1, l2=l2 ,l3=l3)] *
        gamma_hat[matrix_position_to_vector_index_2way(c(j,k), l1=l1, l2=l2 ,l3=l3)] *
        delta_hat[table_position_to_vector_index3(c(i,j,k), l1=l1, l2=l2, l3=l3)]
      
      }}
      
      X1_tilde<- X[,discard_main] + two_ways
      X2_tilde<-three_ways
      
     
      
      
      Q_old <- Q_normal(X=X,y=y, beta=beta_hat, gamma_vec=gamma_hat, delta_vec=delta_hat, 
                        lambda_beta=lambda_beta, lambda_gamma=0, lambda_delta=0, 
                        w_beta=1, w_gamma=1, w_delta=1,l1=l1,l2=l2,l3=l3, already_multiplied=TRUE)
      
      coefs<-get_coef_from_xyz(x=array(X1_tilde), y= array(y_tilde), z=array(X2_tilde) ) #coefs as c0 c1 c2...c4
      beta_hat[k]<-poly_lasso_min(coefs = coefs, lambda = lambda_beta) #beta updated
      #update beta_23way
      beta_2way <- get_beta_vec_2way(beta = beta_hat, l1 = l1, l2 = l2, l3 = l3, gamma = gamma_hat, only_beta = FALSE) ###This is with delta
      beta_3way <- get_beta_vec_3way(beta_2way = beta_2way, l1 = l1, l2 = l2, l3 = l3, delta = delta_hat, only_beta = FALSE) #This is with gamma WITH delta
      
      
      
      Q_new <- Q_normal(X=X,y=y, beta=beta_hat, gamma_vec=gamma_hat, delta_vec=delta_hat, 
                        lambda_beta=lambda_beta, lambda_gamma=0, lambda_delta=0, 
                        w_beta=1, w_gamma=1, w_delta=1,l1=l1,l2=l2,l3=l3, already_multiplied=TRUE)
      if (Q_new-Q_old >= Q_old/100)
      {print("There might be numerical instability in beta.")}
      
      #cat(" new-old: ",Q_new-Q_old, " Q: ",Q_new)
    }

  #print(beta_hat)
  return(beta_hat)
  
  }
  
  
  



data<- create_basic_dataset()
X<- data$X
y<- data$y$obs

beta_true<- data$beta[-1,]
l1=8
l2=8
l3=4
#print(beta_true)
beta_main<-beta_true[1:(l1+l2+l3)]
beta_2way<-beta_true[unlist(get_ranges(l1,l2,l3)[2])]
beta_3way<-beta_true[unlist(get_ranges(l1,l2,l3)[3])]
beta_2way_without_gamma<-get_beta_vec_2way(beta = beta_main, l1=l1, l2=l2, l3=l3, gamma=NULL, only_beta = TRUE )
beta_3way_without_gamma<-get_beta_vec_3way(beta_2way = beta_2way, l1=l1, l2=l2, l3=l3, delta = NULL, only_beta = TRUE)


gamma_true<-beta_2way/beta_2way_without_gamma
gamma_true[is.nan(gamma_true)]<-0
delta_true<-beta_3way/beta_3way_without_gamma
delta_true[is.nan(delta_true)]<-0


delta_hat<-delta_true
gamma_hat<-gamma_true

beta_hat<-beta_main+rnorm(length(beta_main), 0,0.5)
beta_hat[1]<-30

beta_pred <- update_beta(X=X, y=y,beta_hat=beta_hat, gamma_hat=gamma_hat, delta_hat=delta_hat, lambda_beta=100,  l1=l1, l2=l2, l3=l3, w=1) 
beta_hat[17:20]
beta_pred[17:20]
beta_main[17:20]

sum(beta_pred==0)  
sum(abs(beta_pred))


sum(abs(gamma_true))
sum(gamma_true==0)











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
  
  
  
  fit<-function(X, y, lambda_beta, lambda_gamma, lambda_delta, w_beta=NULL, w_gamma=NULL, w_delta=NULL,  tol=1e-5, max_iter=50, compute_Q=Q_normal)
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
    
    
    for (i in c(1:max_iter))  ###print smth IF LOSS DOES NOT DECREASE AFTER ONE ITER
    {    
      ## STEP 2 (UPDATE DELTA)
      delta_hat<- update_delta(X=X, y=y, beta_hat=beta_hat, gamma_hat=gamma_hat, delta_hat=delta_hat, lambda_delta=lambda_delta, l1=self$l1, l2=self$l2, l3=self$l3)
      
      ## STEP 3 (UPDATE GAMMA)
      gamma_hat<- update_gamma(X=X, y=y,beta_hat=beta_hat, gamma_hat=gamma_hat, delta_hat=delta_hat, lambda_gamma=lambda_gamma,  l1=self$l1, l2=self$l2, l3=self$l3, w=1)

      ## STEP 4 (UPDATE BETA)
      beta_hat<- update_beta(X=X, y=y,beta_hat=beta_hat, gamma_hat=gamma_hat, delta_hat=delta_hat, lambda_beta=lambda_beta,  l1=self$l1, l2=self$l2, l3=self$l3, w=1)

      ## STEP 5 (COMPUTE REL_DIF)
      Q_new<-compute_Q(X=X,y=y, beta= beta_hat, gamma_vec=gamma_hat, delta_vec=delta_hat, 
                       lambda_beta=lambda_beta, lambda_gamma=lambda_gamma, lambda_delta=lambda_delta,
                       w_beta =w_beta, w_gamma=w_gamma, w_delta=w_delta,l1=self$l1, l2=self$l2, l3=self$l3)
      rel_dif<-compute_rel_dif(Q_old=Q_old, Q_new=Q_new)

      
      if (Q_new>Q_old*1.01)
      {print("there is numerical instability overall. ")}
      if(i%%2==1)
      {cat("  Q: ",Q_new  )}
      
      
      
      ## STEP 6 (CHECK IF SMALL ENOUGH TO STOP)
      if (abs(rel_dif)<=tol){
        self$beta_hat<-beta_hat
        self$gamma_hat<-gamma_hat
        self$delta_hat<-delta_hat
        beta_2way <- get_beta_vec_2way(beta = self$beta_hat, l1 = l1, l2 = l2, l3 = l3, gamma = self$gamma_hat, only_beta = FALSE) ###This is with delta
        beta_3way <- get_beta_vec_3way(beta_2way = beta_2way, l1 = l1, l2 = l2, l3 = l3, delta = self$delta_hat, only_beta = FALSE) #This is with gamma WITH delta
        beta_all<-array(c(beta_hat, beta_2way, beta_3way))
        
        return (list("beta_hat" = self$beta_hat, "gamma_hat" = self$gamma_hat , "delta_hat" = self$delta_hat, "beta_all" = beta_all)) }
      Q_old<-Q_new #UPDATE Q_old
      }
      cat("It has not converged. The relative difference between last 2 residuals is:", compute_rel_dif(Q_old=Q_old,Q_new=Q_new) )
      self$beta_hat<-beta_hat
      self$gamma_hat<-gamma_hat
      self$delta_hat<-delta_hat
      beta_2way <- get_beta_vec_2way(beta = self$beta_hat, l1 = l1, l2 = l2, l3 = l3, gamma = self$gamma_hat, only_beta = FALSE) ###This is with delta
      beta_3way <- get_beta_vec_3way(beta_2way = beta_2way, l1 = l1, l2 = l2, l3 = l3, delta = self$delta_hat, only_beta = FALSE) #This is with gamma WITH delta
      beta_all<-array(c(beta_hat, beta_2way, beta_3way))
      
      return (list("beta_hat" = self$beta_hat, "gamma_hat" = self$gamma_hat , "delta_hat" = self$delta_hat, "beta_all" = beta_all)) 
      }

  
  
  
  predict<-function(self, X_new, mean_y,scale=FALSE)
  {if (scale ==TRUE)
    {X<-scale(X_new)}
  beta_2way <- get_beta_vec_2way(beta = self$beta_hat, l1 = l1, l2 = l2, l3 = l3, gamma = self$gamma_hat, only_beta = FALSE) ###This is with delta
  beta_3way <- get_beta_vec_3way(beta_2way = beta_2way, l1 = l1, l2 = l2, l3 = l3, delta = self$delta_hat, only_beta = FALSE) #This is with gamma WITH delta
  beta_all<-array(c(beta_hat, beta_2way, beta_3way))
  y_pred<-  X_new%*%beta_all+mean_y
  return(y_pred)
  }
  
  R2_score<-function(self, X_new, y_true, scale=FALSE, verbose=TRUE)
  {y_pred<-predict(self, X_new, scale=scale, mean_y=mean(y_true))
  if (verbose == TRUE)
  {cat ("r2 score is ", r2(y_true, y_pred))
    cat(length(y_true), ' ', length(y_pred))
    plot(array(y_pred), array(y_true))}
  
  return(r2(y_true, y_pred))}
  
  return(list( fit = fit, predict = predict, R2_score = R2_score, self = self))
    
}
  








#### TEST SHIM CLASS ##########
  
  
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
beta_2way_without_gamma<-get_beta_vec_2way(beta = beta_main, l1=l1, l2=l2, l3=l3, gamma=NULL, only_beta = TRUE )
beta_3way_without_gamma<-get_beta_vec_3way(beta_2way = beta_2way, l1=l1, l2=l2, l3=l3, delta = NULL, only_beta = TRUE)


gamma_true<-beta_2way/beta_2way_without_gamma
gamma_true[is.nan(gamma_true)]<-0
delta_true<-beta_3way/beta_3way_without_gamma
delta_true[is.nan(delta_true)]<-0

y_centered<-y-mean(y)
lasso_model <- glmnet(X, y_centered, alpha = 1, intercept = FALSE, standardize = FALSE, lambda=0.1)
coefs_lasso<-coefficients(lasso_model)[-1]
beta_main_lasso<-coefs_lasso[range_main]
beta_2way_lasso<-coefs_lasso[range_theta]
beta_3way_lasso<-coefs_lasso[range_psi]
beta_2way_lasso_without_gamma<-get_beta_vec_2way(beta_main_lasso,l1=l1,l2=l2,l3=l3,only_beta = TRUE)
beta_3way_lasso_without_delta<- get_beta_vec_3way(beta_2way_lasso, l1=l1, l2=l2, l3=l3, only_beta = TRUE)


length(c(beta_main, beta_2way, beta_3way))
print("lasso")
all_beta_functions(beta_main, beta_main_lasso)
all_beta_functions(beta_2way, beta_2way_lasso)
all_beta_functions(beta_3way, beta_3way_lasso)



beta_hat<-beta_main_lasso
gamma_hat<- beta_2way_lasso/beta_2way_lasso_without_gamma
gamma_hat[is.nan(gamma_hat)]<-0
gamma_hat[!is.finite(gamma_hat)]<-0 #this is 0 in shim case
delta_hat<- beta_3way_lasso/beta_3way_lasso_without_delta
delta_hat[!is.finite(delta_hat)]<-0
delta_hat[is.nan(delta_hat)]<-0

##USE SHIM MODEL

lambda_beta<-100
lambda_gamma<-1200
lambda_delta<-500


my_shim<-SHIM_3way(X=X, y=y, beta_init = beta_hat, gamma_init = gamma_hat, delta_init = delta_hat, l1=l1, l2=l2, l3=l3, scale = FALSE)
fitted<-my_shim$fit(X=X, y=y, lambda_beta = lambda_beta, 
                    lambda_gamma = lambda_gamma, lambda_delta = lambda_delta, w_beta = 1, w_gamma = 1, w_delta = 1, tol=6e-3)

my_shim$R2_score(self=fitted, X_new=X, y_true=y )


beta_all_shim<-fitted$beta_all
beta_main_shim<-beta_all_shim[range_main]
beta_2way_shim<-beta_all_shim[range_theta]
beta_3way_shim<-beta_all_shim[range_psi]



all_beta_functions(beta_main, beta_main_shim)
all_beta_functions(beta_2way, beta_2way_shim)
all_beta_functions(beta_3way, beta_3way_shim)




