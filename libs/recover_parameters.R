
assert <- function(condition, message) {
  if (!condition)
    stop(message)
}

set_0_threshold<- function(x,threshold)
{ if (abs(x)<=threshold)
{return (0)}
  return (x)
}


#Function transform maineff+idx to position theta
get_position_in_theta<-function(main, idx, l1, l2, l3)
{assert(main ==1 | main==2 | main ==3, " main should be 1 2 or 3, i.e. which main effect it is")
  assert(idx<=max(l1,l2,l3), "idx is too big")
  pre_position<-c(0,l1, l1+l2)
  position_theta<-pre_position[main]+ idx
  return(position_theta)
}

get_position_in_theta(3,1,4,4,4)



get_all_beta<- function(beta, l1, l2, l3, threshold=1e-3)
{ assert(length(beta)==l1+l2+l3, "length of coefs is not right")
  beta_old<-beta
 beta_new<-array(0, dim=l1+l2+l3+3) # shape of new beta

 beta_new[1:l1]<-beta_old[1:l1]
 beta_new[l1+1]<- set_0_threshold(x=-sum(beta_old[1:l1]), threshold = threshold)
 
 beta_new[(l1+2):(l1+l2+1)]<-beta_old[(l1+1): (l1+l2)]
 beta_new[l1+l2+2]<- set_0_threshold(x=-sum(beta_old[(l1+1): (l1+l2)]), threshold = threshold) 
 
 beta_new[(l1+l2+3):(l1+l2+l3+2)]<-beta_old[(l1+l2+1): (l1+l2+l3)]
 beta_new[l1+l2+l3+3]<- set_0_threshold(x = -sum(beta_old[(l1+l2+1): (l1+l2+l3)]), threshold = threshold)
 
 return(beta_new)
}

#get_all_beta(beta=c(1,2,3,4,5,6,7),l1=3, l2=2, l3=2, threshold = 1e1)

get_all_theta <-function(theta, l1, l2, l3, threshold=1e-3)
{ theta_old<-theta

range1_old<-c(1:l1)
range2_old<-c( (l1+1): (l1+l2)  )
range3_old<-c( (l1+l2+1) : (l1+l2+l3)   )
range1_new<-c(1:(l1+1))
range2_new<-c( (l1+2): (l1+l2+2)  )
range3_new<-c( (l1+l2+3) : (l1+l2+l3+3)   )


ls_range_old<-list(range1_old, range2_old, range3_old)
ls_range_new<-list(range1_new, range2_new, range3_new)

#Step1 create theta new
theta_new<- matrix(0, nrow = l1+l2+l3+3, ncol=l1+l2+l3+3) 

for (i in c(1:3))
{for (j in c(1:3))
{rangenew_i<-unlist(ls_range_new[i])
 rangenew_j<-unlist(ls_range_new[j])
 rangeold_i<-unlist(ls_range_old[i])
 rangeold_j<-unlist(ls_range_old[j])
    theta_new[rangenew_i[-length(rangenew_i)], rangenew_j[-length(rangenew_j)]] <- theta_old[rangeold_i, rangeold_j]}}

print(theta_new)
##STEP 2 Create matrix of coefs
theta<- (theta+t(theta) )/2

#STEP 3 Identify (ab)_{il1} ...


#STEP 4 IDENTIFY (ab)_{l1l2}...


  
return(theta_new)  
  
}

theta<-matrix(1, nrow=7,ncol=7)
get_all_theta(theta, l1=2,l2=2,l3=3 )


p