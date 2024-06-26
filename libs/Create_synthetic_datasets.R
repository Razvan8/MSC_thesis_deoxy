# Clear all
#rm(list=ls())
#gc()
#try(dev.off(dev.list()["RStudioGD"]),silent=TRUE)
#try(dev.off(),silent=TRUE)
#cat("\014")

# Memory
#memory.limit(64000)

set.seed(123)

#######################
## Dummy Function
#######################

dummy.matrix <- function(NF=2, NL=rep(2,NF)) {
  # Computes dummy matrix from number of factors NF and number of levels NL
  # NF is an integer between 2 and 4
  # NL is a vector of length NF having integer entries between 2 and 100
  
  # Factors and Levels
  fac.tags <- c("A", "B", "C", "D")
  fac.levs <- as.character(1:100)
  
  if (NF == 2) {
    
    # One-way
    L.1 <- NL[1]
    L.2 <- NL[2]
    fac.1 <- paste(fac.tags[1], fac.levs[1:L.1], sep=".")
    fac.2 <- paste(fac.tags[2], fac.levs[1:L.2], sep=".")
    
    # Two-ways
    L.12 <- L.1*L.2
    fac.12 <- sort(as.vector(outer(fac.1, fac.2, paste, sep=":")))
    
    # Dummy design matrix 2-way
    n.2w <- L.12
    p.2w <- n.2w-1
    x.2w <- data.frame(matrix(0, nrow=n.2w, ncol=p.2w))
    rownames(x.2w) <- c(fac.12)
    colnames(x.2w) <- c(fac.1[-L.1], fac.2[-L.2],
                        sort(as.vector(outer(fac.1[-L.1], fac.2[-L.2], paste, sep=":")))
    )
    for (col in 1:ncol(x.2w)) {
      col.tags <- unlist(strsplit(colnames(x.2w)[col], split=":"))
      if (length(col.tags)==1) {
        fac.tag <- unlist(strsplit(col.tags,split="\\."))
        x.2w[grepl(col.tags,split(unlist(strsplit(rownames(x.2w),split=":")), 1:2)[[1]]),col] <- 1
        x.2w[grepl(col.tags,split(unlist(strsplit(rownames(x.2w),split=":")), 1:2)[[2]]),col] <- 1
        x.2w[grepl(paste(fac.tag[1],L.1,sep="."),split(unlist(strsplit(rownames(x.2w),split=":")), 1:2)[[1]]),col] <- -1
        x.2w[grepl(paste(fac.tag[1],L.2,sep="."),split(unlist(strsplit(rownames(x.2w),split=":")), 1:2)[[2]]),col] <- -1
      }
      if (length(col.tags)==2) {
        col.1 <- which(colnames(x.2w)==col.tags[1])
        col.2 <- which(colnames(x.2w)==col.tags[2])
        x.2w[,col] <- x.2w[,col.1]*x.2w[,col.2] 
      }
    }
    
    return(x.2w)
  }
  
  if (NF == 3) {
    
    # One-way
    L.1 <- NL[1]
    L.2 <- NL[2]
    L.3 <- NL[3]
    fac.1 <- paste(fac.tags[1], fac.levs[1:L.1], sep=".")
    fac.2 <- paste(fac.tags[2], fac.levs[1:L.2], sep=".")
    fac.3 <- paste(fac.tags[3], fac.levs[1:L.3], sep=".")
    
    # Two-ways
    L.12 <- L.1*L.2
    L.13 <- L.1*L.3
    L.23 <- L.2*L.3
    fac.12 <- sort(as.vector(outer(fac.1, fac.2, paste, sep=":")))
    fac.13 <- sort(as.vector(outer(fac.1, fac.3, paste, sep=":")))
    fac.23 <- sort(as.vector(outer(fac.2, fac.3, paste, sep=":")))
    
    # Three-ways
    L.123 <- L.1*L.2*L.3
    fac.123 <- sort(as.vector(outer(fac.1, fac.23, paste, sep=":")))
    
    # Dummy design matrix 3-way
    n.3w <- L.123
    p.3w <- n.3w-1
    x.3w <- data.frame(matrix(0, nrow=n.3w, ncol=p.3w))
    rownames(x.3w) <- c(fac.123)
    colnames(x.3w) <- c(fac.1[-L.1], fac.2[-L.2], fac.3[-L.3],
                        sort(as.vector(outer(fac.1[-L.1], fac.2[-L.2], paste, sep=":"))),
                        sort(as.vector(outer(fac.1[-L.1], fac.3[-L.3], paste, sep=":"))),
                        sort(as.vector(outer(fac.2[-L.2], fac.3[-L.3], paste, sep=":"))),
                        sort(as.vector(outer(fac.1[-L.1], sort(as.vector(outer(fac.2[-L.2], fac.3[-L.3], paste, sep=":"))), paste, sep=":")))
    )
    for (col in 1:ncol(x.3w)) {
      col.tags <- unlist(strsplit(colnames(x.3w)[col], split=":"))
      if (length(col.tags)==1) {
        fac.tag <- unlist(strsplit(col.tags,split="\\."))
        x.3w[grepl(col.tags,split(unlist(strsplit(rownames(x.3w),split=":")), 1:3)[[1]]),col] <- 1
        x.3w[grepl(col.tags,split(unlist(strsplit(rownames(x.3w),split=":")), 1:3)[[2]]),col] <- 1
        x.3w[grepl(col.tags,split(unlist(strsplit(rownames(x.3w),split=":")), 1:3)[[3]]),col] <- 1
        x.3w[grepl(paste(fac.tag[1],L.1,sep="."),split(unlist(strsplit(rownames(x.3w),split=":")), 1:3)[[1]]),col] <- -1
        x.3w[grepl(paste(fac.tag[1],L.2,sep="."),split(unlist(strsplit(rownames(x.3w),split=":")), 1:3)[[2]]),col] <- -1
        x.3w[grepl(paste(fac.tag[1],L.3,sep="."),split(unlist(strsplit(rownames(x.3w),split=":")), 1:3)[[3]]),col] <- -1
      }
      if (length(col.tags)==2) {
        col.1 <- which(colnames(x.3w)==col.tags[1])
        col.2 <- which(colnames(x.3w)==col.tags[2])
        x.3w[,col] <- x.3w[,col.1]*x.3w[,col.2]
      }
      if (length(col.tags)==3) {
        col.1 <- which(colnames(x.3w)==col.tags[1])
        col.2 <- which(colnames(x.3w)==col.tags[2])
        col.3 <- which(colnames(x.3w)==col.tags[3])
        x.3w[,col] <- x.3w[,col.1]*x.3w[,col.2]*x.3w[,col.3]
      }
    }
    
    return(x.3w)
  }
  
  if (NF == 4) {
    
    # One-way
    L.1 <- NL[1]
    L.2 <- NL[2]
    L.3 <- NL[3]
    L.4 <- NL[4]
    fac.1 <- paste(fac.tags[1], fac.levs[1:L.1], sep=".")
    fac.2 <- paste(fac.tags[2], fac.levs[1:L.2], sep=".")
    fac.3 <- paste(fac.tags[3], fac.levs[1:L.3], sep=".")
    fac.4 <- paste(fac.tags[4], fac.levs[1:L.4], sep=".")
    
    # Two-ways
    L.12 <- L.1*L.2
    L.13 <- L.1*L.3
    L.14 <- L.1*L.4
    L.23 <- L.2*L.3
    L.24 <- L.2*L.4
    L.34 <- L.3*L.4
    fac.12 <- sort(as.vector(outer(fac.1, fac.2, paste, sep=":")))
    fac.13 <- sort(as.vector(outer(fac.1, fac.3, paste, sep=":")))
    fac.14 <- sort(as.vector(outer(fac.1, fac.4, paste, sep=":")))
    fac.23 <- sort(as.vector(outer(fac.2, fac.3, paste, sep=":")))
    fac.24 <- sort(as.vector(outer(fac.2, fac.4, paste, sep=":")))
    fac.34 <- sort(as.vector(outer(fac.3, fac.4, paste, sep=":")))
    
    # Three-ways
    L.123 <- L.1*L.2*L.3
    L.124 <- L.1*L.2*L.4
    L.134 <- L.1*L.3*L.4
    L.234 <- L.2*L.3*L.4
    fac.123 <- sort(as.vector(outer(fac.1, fac.23, paste, sep=":")))
    fac.124 <- sort(as.vector(outer(fac.1, fac.24, paste, sep=":")))
    fac.134 <- sort(as.vector(outer(fac.1, fac.34, paste, sep=":")))
    fac.234 <- sort(as.vector(outer(fac.2, fac.34, paste, sep=":")))
    
    # Four-ways
    L.1234 <- L.1*L.2*L.3*L.4
    fac.1234 <- sort(as.vector(outer(fac.1, fac.234, paste, sep=":")))
    
    # Dummy design matrix 4-way
    n.4w <- L.1234
    p.4w <- n.4w-1
    x.4w <- data.frame(matrix(0, nrow=n.4w, ncol=p.4w))
    rownames(x.4w) <- c(fac.1234)
    colnames(x.4w) <- c(fac.1[-L.1], fac.2[-L.2], fac.3[-L.3], fac.4[-L.4],
                        sort(as.vector(outer(fac.1[-L.1], fac.2[-L.2], paste, sep=":"))),
                        sort(as.vector(outer(fac.1[-L.1], fac.3[-L.3], paste, sep=":"))),
                        sort(as.vector(outer(fac.1[-L.1], fac.4[-L.4], paste, sep=":"))),
                        sort(as.vector(outer(fac.2[-L.2], fac.3[-L.3], paste, sep=":"))),
                        sort(as.vector(outer(fac.2[-L.2], fac.4[-L.4], paste, sep=":"))),
                        sort(as.vector(outer(fac.3[-L.3], fac.4[-L.4], paste, sep=":"))),
                        sort(as.vector(outer(fac.1[-L.1], sort(as.vector(outer(fac.2[-L.2], fac.3[-L.3], paste, sep=":"))), paste, sep=":"))),
                        sort(as.vector(outer(fac.1[-L.1], sort(as.vector(outer(fac.2[-L.2], fac.4[-L.4], paste, sep=":"))), paste, sep=":"))),
                        sort(as.vector(outer(fac.1[-L.1], sort(as.vector(outer(fac.3[-L.3], fac.4[-L.4], paste, sep=":"))), paste, sep=":"))),
                        sort(as.vector(outer(fac.2[-L.2], sort(as.vector(outer(fac.3[-L.3], fac.4[-L.4], paste, sep=":"))), paste, sep=":"))),
                        sort(as.vector(outer(fac.1[-L.1], sort(as.vector(outer(fac.2[-L.2], sort(as.vector(outer(fac.3[-L.3], fac.4[-L.4], paste, sep=":"))),
                                                                               paste, sep=":"))), paste, sep=":")))
    )
    for (col in 1:ncol(x.4w)) {
      col.tags <- unlist(strsplit(colnames(x.4w)[col], split=":"))
      if (length(col.tags)==1) {
        fac.tag <- unlist(strsplit(col.tags,split="\\."))
        x.4w[grepl(col.tags,split(unlist(strsplit(rownames(x.4w),split=":")), 1:4)[[1]]),col] <- 1
        x.4w[grepl(col.tags,split(unlist(strsplit(rownames(x.4w),split=":")), 1:4)[[2]]),col] <- 1
        x.4w[grepl(col.tags,split(unlist(strsplit(rownames(x.4w),split=":")), 1:4)[[3]]),col] <- 1
        x.4w[grepl(col.tags,split(unlist(strsplit(rownames(x.4w),split=":")), 1:4)[[4]]),col] <- 1
        x.4w[grepl(paste(fac.tag[1],L.1,sep="."),split(unlist(strsplit(rownames(x.4w),split=":")), 1:4)[[1]]),col] <- -1
        x.4w[grepl(paste(fac.tag[1],L.2,sep="."),split(unlist(strsplit(rownames(x.4w),split=":")), 1:4)[[2]]),col] <- -1
        x.4w[grepl(paste(fac.tag[1],L.3,sep="."),split(unlist(strsplit(rownames(x.4w),split=":")), 1:4)[[3]]),col] <- -1
        x.4w[grepl(paste(fac.tag[1],L.4,sep="."),split(unlist(strsplit(rownames(x.4w),split=":")), 1:4)[[4]]),col] <- -1
      }
      if (length(col.tags)==2) {
        col.1 <- which(colnames(x.4w)==col.tags[1])
        col.2 <- which(colnames(x.4w)==col.tags[2])
        x.4w[,col] <- x.4w[,col.1]*x.4w[,col.2]
      }
      if (length(col.tags)==3) {
        col.1 <- which(colnames(x.4w)==col.tags[1])
        col.2 <- which(colnames(x.4w)==col.tags[2])
        col.3 <- which(colnames(x.4w)==col.tags[3])
        x.4w[,col] <- x.4w[,col.1]*x.4w[,col.2]*x.4w[,col.3]
      }
      if (length(col.tags)==4) {
        col.1 <- which(colnames(x.4w)==col.tags[1])
        col.2 <- which(colnames(x.4w)==col.tags[2])
        col.3 <- which(colnames(x.4w)==col.tags[3])
        col.4 <- which(colnames(x.4w)==col.tags[4])
        x.4w[,col] <- x.4w[,col.1]*x.4w[,col.2]*x.4w[,col.3]*x.4w[,col.4]
      }
    }
    
    return(x.4w)
  }
  
}



#######################
## Toy Model
#######################


libs_path<-file.path("..","libs")
source(file.path(libs_path,'helper_functions.R'))


## Mini hierachical dataset



## Small sparse hierarchical dataset
create_sparse_hier_dataset<-function()
{

x.3w <- dummy.matrix(NF=3, NL=c(5,10,4))

beta<-array(0,dim(x.3w)[2])

#main effects
beta[1]=1 #A1
beta[4]=4 #A4
beta[9]=9 #B5
beta[13]=-13 #B9
beta[16]=-1 #C3

#2 way
beta[21]=2 #A1:B5
beta[25]=-2 #A1:B9
beta[91]=9 #B9:C3

#3way
beta[106]=2  #A1B5C3
beta[105]=1  #A1B5C2
beta[118]=-3  #A1B9C3

noise<-rnorm(dim(x.3w)[1],0, 0.1 )
X<-as.matrix(x.3w)
y<-X%*%beta+noise


# Function to generate Y with a specific SNR
generate_y_with_snr <- function(X, beta, snr) {
  # Calculate the signal

  signal <- X %*% beta
  # Calculate the variance of the signal
  signal_variance <- var(as.vector(signal))
  # Calculate the variance of the noise
  noise_variance <- signal_variance / snr
  # Generate the noise
  epsilon <- rnorm(nrow(X), mean = 0, sd = sqrt(noise_variance))
  print("noise var")
  print(noise_variance)
  # Generate the response
  Y <- signal + epsilon
  return(Y)
}

# Generate data with different SNR values
Y_high_snr <- generate_y_with_snr(X, beta, snr = 10)  # High SNR
Y_moderate_snr <- generate_y_with_snr(X, beta, snr = 1.5)  # Moderate SNR
Y_low_snr <- generate_y_with_snr(X, beta, snr = 1)  # Low SNR

# Check the variances and SNRs
signal_variance <- var(as.vector(X %*% beta))
high_snr <- signal_variance / var(Y_high_snr - X %*% beta)
moderate_snr <- signal_variance / var(Y_moderate_snr - X %*% beta)
low_snr <- signal_variance / var(Y_low_snr - X %*% beta)

print(paste("High SNR:", high_snr))
print(paste("Moderate SNR:", moderate_snr))
print(paste("Low SNR:", low_snr))




return(list('X'=X, 'y'=y, 'beta'=beta))

}

#data=create_sparse_hier_dataset()
#X<-data$X
#y<-data$y
#print(dim(X))
#print(dim(y))




##Hierarchical dataset ############################## quite sparse##

create_hier_dataset_medium<-function(){
x.3w <- dummy.matrix(NF=3, NL=c(6,5,4))
# Hierarchical Coefficients (2way)
p.3w <- ncol(x.3w)
n.3w <- p.3w + 1
beta.min <- 1
beta.max <- 10
beta.true <- data.frame(rep(0, n.3w))
rownames(beta.true) <- c("interc", colnames(x.3w))
colnames(beta.true) <- c("coeffs")
beta.true$coeffs <- runif(n.3w, beta.min, beta.max)*sample(c(1,-1),size=n.3w,replace=TRUE)

levs.true <- c("interc","A.1", "A.2", "A.3","A.4",  "B.1", "B.2","B.3", "C.1", "C.2","C.3",  "A.1:B.1","A.1:B.2",
               "A.2:B.1","A.2:B.2", "A.3:B.1","A.3:B.2",
               "A.1:C.1","A.1:C.2","A.2:C.1", "A.2:C.2","A.3:C.1","A.3:C.2", 
               "B.1:C1","B.1:C2", "B2:C1", "B.2:C2", 
               "A.1:B.1:C.1","A.1:B.1:C.2","A.1:B.2:C.1", "A.1:B.2:C.2","A.2:B.1:C.1","A.2:B.1:C.2","A.2:B.2:C.2","A.3:B.1:C.1" )

beta.true$coeffs[which(!is.element(rownames(beta.true),levs.true))] <- 0
#beta.true

# Response vector (2way)
sigma.y <- 3
y.3w <- data.frame(row.names=rownames(x.3w))
y.3w$obs <- beta.true$coeffs[1] + as.matrix(x.3w)%*%as.vector(beta.true$coeffs)[-1] + rnorm(nrow(y.3w), 0, sigma.y)
y.3w$true <- beta.true$coeffs[1] + as.matrix(x.3w)%*%as.vector(beta.true$coeffs)[-1]
return(list('X'=as.matrix(x.3w), 'y'=y.3w, 'beta'=beta.true))
}



create_hier_dataset_medium_2way<-function(){
  x.3w <- dummy.matrix(NF=3, NL=c(6,5,4))

  p.3w <- ncol(x.3w)
  n.3w <- p.3w + 1
  beta.min <- 1
  beta.max <- 10
  beta.true <- data.frame(rep(0, n.3w))
  rownames(beta.true) <- c("interc", colnames(x.3w))
  colnames(beta.true) <- c("coeffs")
  beta.true$coeffs <- runif(n.3w, beta.min, beta.max)*sample(c(1,-1),size=n.3w,replace=TRUE)
  
  levs.true <- c("interc","A.1", "A.2", "A.3","A.4",  "B.1", "B.2","B.3", "C.1", "C.2","C.3",  "A.1:B.1","A.1:B.2",
                 "A.2:B.1","A.2:B.2", "A.3:B.1","A.3:B.2",
                 "A.1:C.1","A.1:C.2","A.2:C.1", "A.2:C.2","A.3:C.1","A.3:C.2", 
                 "B.1:C1","B.1:C2", "B2:C1", "B.2:C2")
  
  beta.true$coeffs[which(!is.element(rownames(beta.true),levs.true))] <- 0
  #beta.true
  
  get_sigma_from_snr <- function(X, beta, snr) {
    # Calculate the signal
    
    signal <- X %*% beta
    # Calculate the variance of the signal
    signal_variance <- var(as.vector(signal))
    # Calculate the variance of the noise
    noise_variance <- signal_variance / snr
    cat("Error sigma: ", srqt(noise_variance))
    print(" ")
    return(sqrt(noise_variance))
  }
  
  # Generate data with different SNR values
  Y_high_snr <- generate_y_with_snr(X, beta, snr = 10)  # High SNR
  Y_moderate_snr <- generate_y_with_snr(X, beta, snr = 1.5)  # Moderate SNR
  Y_low_snr <- generate_y_with_snr(X, beta, snr = 1)  # Low SNR
  
  # Check the variances and SNRs
  signal_variance <- var(as.vector(X %*% beta))
  high_snr <- signal_variance / var(Y_high_snr - X %*% beta)
  moderate_snr <- signal_variance / var(Y_moderate_snr - X %*% beta)
  low_snr <- signal_variance / var(Y_low_snr - X %*% beta)
  
  print(paste("High SNR:", high_snr))
  print(paste("Moderate SNR:", moderate_snr))
  print(paste("Low SNR:", low_snr))
  
  
  sigma.y <- get_sigma_from_snr(x.3w, beta.true,snr=2)
  y.3w <- data.frame(row.names=rownames(x.3w))
  y.3w$obs <- beta.true$coeffs[1] + as.matrix(x.3w)%*%as.vector(beta.true$coeffs)[-1] + rnorm(nrow(y.3w), 0, sigma.y)
  y.3w$true <- beta.true$coeffs[1] + as.matrix(x.3w)%*%as.vector(beta.true$coeffs)[-1]
  return(list('X'=as.matrix(x.3w), 'y'=y.3w, 'beta'=beta.true))
}








##Hierarchical dataset similar to paper ############################################################################


create_hier_dataset_paper<-function(){
  x.3w <- dummy.matrix(NF=3, NL=c(15,10,5))
  # Hierarchical Coefficients (2way)
  p.3w <- ncol(x.3w)
  n.3w <- p.3w + 1
  beta.min <- 1
  beta.max <- 10
  beta.true <- data.frame(rep(0, n.3w))
  rownames(beta.true) <- c("interc", colnames(x.3w))
  colnames(beta.true) <- c("coeffs")
  beta.true$coeffs <- runif(n.3w, beta.min, beta.max)*sample(c(1,-1),size=n.3w,replace=TRUE)
  
  levs.true <- c("interc","A.1", "A.2", "A.3","A.4",  "B.1", "B.2","B.3", "C.1", "C.2","C.3",
                 "A.1:B.1","A.1:B.2","A.2:B.1","A.2:B.2",
                 "A.1:C.1","A.1:C.2", "A.2:C.1","A.2:C.2",
                 "B.1:C.1","B.3:C.3",
                 "A.1:B.1:C.1","A.1:B.1:C.2","A.1:B.2:C.1", "A.1:B.2:C.2","A.2:B.1:C.1","A.2:B.1:C.2","A.2:B.2:C.2" )
  
  beta.true$coeffs[which(!is.element(rownames(beta.true),levs.true))] <- 0
  #beta.true
  
  # Response vector (2way)
  sigma.y <- 3
  y.3w <- data.frame(row.names=rownames(x.3w))
  y.3w$obs <- beta.true$coeffs[1] + as.matrix(x.3w)%*%as.vector(beta.true$coeffs)[-1] + rnorm(nrow(y.3w), 0, sigma.y)
  y.3w$true <- beta.true$coeffs[1] + as.matrix(x.3w)%*%as.vector(beta.true$coeffs)[-1]
  return(list('X'=as.matrix(x.3w), 'y'=y.3w, 'beta'=beta.true))
}



set.seed(123)
create_hier_dataset_paper_3way<-function(){
  set.seed(123)
  x.3w <- dummy.matrix(NF=3, NL=c(10,10,5))
  l1=9
  l2=9
  l3=4
  col_mains<-colnames(x.3w)[c(1:(l1+l2+l3))]
  #print(col_mains)
  col_psi<-colnames(x.3w)[c ( (l1+l2+l3+ l1*l2+ l2*l3 + l1*l3 + 1): (l1+l2+l3+ l1*l2+ l2*l3+ +l1*l3+ l1*l2*l3) )]
  col_theta_good<-c()
   for (i in c(1:l1)) {
      for (j in c(1:l2)) {
        # Create the string and append to the vector
        col_theta_good <- c(col_theta_good, paste0("A.", i, ":B.", j))
      }
     for(k in c(1:l3))
     {col_theta_good <- c(col_theta_good, paste0("A.", i, ":C.", k))}
   }
  for (j in c(1:l2))
  {for(k in c(1:l3))
  {col_theta_good <- c(col_theta_good, paste0("B.", j, ":C.", k))}
  }
  #print(col_theta_good)
  col_all_good=c(col_mains,col_theta_good,col_psi)
  #print(col_all_good)
  # Hierarchical Coefficients (2way)
  p.3w <- ncol(x.3w)
  n.3w <- p.3w + 1
  beta.min <- 1
  beta.max <- 10
  beta.true <- data.frame(rep(0, n.3w))
  rownames(beta.true) <- c("interc", colnames(x.3w))
  colnames(beta.true) <- c("coeffs")
  set.seed(123)
  beta.true$coeffs <- runif(n.3w, beta.min, beta.max)*sample(c(1,-1),size=n.3w,replace=TRUE)
  
  
  levs.true <- c("interc","A.1", "A.2", "A.3","A.4",  "B.1", "B.2","B.3", "C.1", "C.2","C.3",
                 "A.1:B.1","A.1:B.2","A.2:B.1","A.2:B.2",
                 "A.1:C.1","A.1:C.2", "A.2:C.1","A.2:C.2",
                 "B.1:C.1","B.3:C.3",
                 "A.1:B.1:C.1","A.1:B.1:C.2","A.1:B.2:C.1", "A.1:B.2:C.2","A.2:B.1:C.1","A.2:B.1:C.2")
  
  beta.true$coeffs[which(!is.element(rownames(beta.true),levs.true))] <- 0
  beta.true_new <- beta.true[c("interc", col_all_good), , drop=FALSE ]
  #print("new beta true")
  #print(beta.true_new)
  rownames(beta.true_new)<-c("interc", col_all_good)
  #beta.true
  #print(beta.true)
  # Response vector (2way)
  sigma.y <- 5
  y.3w <- data.frame(row.names=rownames(x.3w))
  set.seed(123)
  y.3w$obs <- beta.true$coeffs[1] + as.matrix(x.3w)%*%as.vector(beta.true$coeffs)[-1] + rnorm(nrow(y.3w), 0, sigma.y)
  y.3w$true <- beta.true$coeffs[1] + as.matrix(x.3w)%*%as.vector(beta.true$coeffs)[-1]
  #print(colnames(x.3w))
  x.3w_new<-x.3w[,col_all_good ,drop=FALSE]
  
  if(all(rownames(beta.true_new)[-1]==colnames(x.3w_new)) ==TRUE)
  {print("Data loaded properly")}
  print("SNR:")
  print(var(as.vector(y.3w$true))/ (sigma.y^2) )
  #print("newcolnames")
  #print(colnames(x.3w_new))
  return(list('X'=as.matrix(x.3w_new), 'y'=y.3w, 'beta'=beta.true_new))
}




set.seed(123)
create_hier_dataset_paper_many_main<-function(){
  set.seed(123)
  x.3w <- dummy.matrix(NF=3, NL=c(7,6,5))
  l1=6
  l2=5
  l3=4
  col_mains<-colnames(x.3w)[c(1:(l1+l2+l3))]
  #print(col_mains)
  col_psi<-colnames(x.3w)[c ( (l1+l2+l3+ l1*l2+ l2*l3 + l1*l3 + 1): (l1+l2+l3+ l1*l2+ l2*l3+ +l1*l3+ l1*l2*l3) )]
  col_theta_good<-c()
  for (i in c(1:l1)) {
    for (j in c(1:l2)) {
      # Create the string and append to the vector
      col_theta_good <- c(col_theta_good, paste0("A.", i, ":B.", j))
    }
    for(k in c(1:l3))
    {col_theta_good <- c(col_theta_good, paste0("A.", i, ":C.", k))}
  }
  for (j in c(1:l2))
  {for(k in c(1:l3))
  {col_theta_good <- c(col_theta_good, paste0("B.", j, ":C.", k))}
  }
  #print(col_theta_good)
  col_all_good=c(col_mains,col_theta_good,col_psi)
  #print(col_all_good)
  # Hierarchical Coefficients (2way)
  p.3w <- ncol(x.3w)
  n.3w <- p.3w + 1
  beta.min <- 1
  beta.max <- 10
  beta.true <- data.frame(rep(0, n.3w))
  rownames(beta.true) <- c("interc", colnames(x.3w))
  colnames(beta.true) <- c("coeffs")
  set.seed(123)
  beta.true$coeffs <- runif(n.3w, beta.min, beta.max)*sample(c(1,-1),size=n.3w,replace=TRUE)
  
  
  levs.true <- c("interc","A.1", "A.2", "A.3","A.4",  "B.1", "B.2","B.3", "C.1", "C.2","C.3",
                 "A.1:B.1","A.1:B.2","A.2:B.1","A.2:B.2",
                 "A.1:C.1","A.1:C.2", "A.2:C.1","A.2:C.2",
                 "B.1:C.1","B.3:C.3",
                 "A.1:B.1:C.1","A.1:B.1:C.2","A.1:B.2:C.1", "A.1:B.2:C.2","A.2:B.1:C.1","A.2:B.1:C.2")
  
  beta.true$coeffs[which(!is.element(rownames(beta.true),levs.true))] <- 0
  beta.true_new <- beta.true[c("interc", col_all_good), , drop=FALSE ]
  #print("new beta true")
  #print(beta.true_new)
  rownames(beta.true_new)<-c("interc", col_all_good)
  #beta.true
  #print(beta.true)
  # Response vector (2way)
  sigma.y <- 5
  y.3w <- data.frame(row.names=rownames(x.3w))
  set.seed(123)
  y.3w$obs <- beta.true$coeffs[1] + as.matrix(x.3w)%*%as.vector(beta.true$coeffs)[-1] + rnorm(nrow(y.3w), 0, sigma.y)
  y.3w$true <- beta.true$coeffs[1] + as.matrix(x.3w)%*%as.vector(beta.true$coeffs)[-1]
  #print(colnames(x.3w))
  x.3w_new<-x.3w[,col_all_good ,drop=FALSE]
  
  if(all(rownames(beta.true_new)[-1]==colnames(x.3w_new)) ==TRUE)
  {print("Data loaded properly")}
  print("SNR:")
  print(var(as.vector(y.3w$true))/ (sigma.y^2) )
  #print("newcolnames")
  #print(colnames(x.3w_new))
  return(list('X'=as.matrix(x.3w_new), 'y'=y.3w, 'beta'=beta.true_new))
}






set.seed(123)
create_hier_dataset_paper_many_main_vary_interactions<-function(magnitude_scale=3){
  set.seed(123)
  x.3w <- dummy.matrix(NF=3, NL=c(7,6,5))
  l1=6
  l2=5
  l3=4
  col_mains<-colnames(x.3w)[c(1:(l1+l2+l3))]
  #print(col_mains)
  col_psi<-colnames(x.3w)[c ( (l1+l2+l3+ l1*l2+ l2*l3 + l1*l3 + 1): (l1+l2+l3+ l1*l2+ l2*l3+ +l1*l3+ l1*l2*l3) )]
  col_theta_good<-c()
  for (i in c(1:l1)) {
    for (j in c(1:l2)) {
      # Create the string and append to the vector
      col_theta_good <- c(col_theta_good, paste0("A.", i, ":B.", j))
    }
    for(k in c(1:l3))
    {col_theta_good <- c(col_theta_good, paste0("A.", i, ":C.", k))}
  }
  for (j in c(1:l2))
  {for(k in c(1:l3))
  {col_theta_good <- c(col_theta_good, paste0("B.", j, ":C.", k))}
  }
  #print(col_theta_good)
  col_all_good=c(col_mains,col_theta_good,col_psi)
  #print(col_all_good)
  # Hierarchical Coefficients (2way)
  p.3w <- ncol(x.3w)
  n.3w <- p.3w + 1
  beta.min <- 1
  beta.max <- 10
  beta.true <- data.frame(rep(0, n.3w))
  rownames(beta.true) <- c("interc", colnames(x.3w))
  colnames(beta.true) <- c("coeffs")
  set.seed(123)
  beta.true$coeffs <- runif(n.3w, beta.min, beta.max)*sample(c(1,-1),size=n.3w,replace=TRUE)
  
  
  levs.true <- c("interc","A.1", "A.2", "A.3","A.4",  "B.1", "B.2","B.3", "C.1", "C.2","C.3",
                 "A.1:B.1","A.1:B.2","A.2:B.1","A.2:B.2",
                 "A.1:C.1","A.1:C.2", "A.2:C.1","A.2:C.2",
                 "B.1:C.1","B.3:C.3",
                 "A.1:B.1:C.1","A.1:B.1:C.2","A.1:B.2:C.1", "A.1:B.2:C.2","A.2:B.1:C.1","A.2:B.1:C.2")
  levs.high<-c("A.1:B.1","A.1:B.2","A.2:B.1","A.2:B.2",
               "A.1:C.1","A.1:C.2", "A.2:C.1","A.2:C.2",
               "B.1:C.1","B.3:C.3",
               "A.1:B.1:C.1","A.1:B.1:C.2","A.1:B.2:C.1", "A.1:B.2:C.2","A.2:B.1:C.1","A.2:B.1:C.2")
  
  beta.true$coeffs[which(!is.element(rownames(beta.true),levs.true))] <- 0
  beta.true$coeffs[which(is.element(rownames(beta.true),levs.high))]<- magnitude_scale*beta.true$coeffs[which(is.element(rownames(beta.true),levs.high))]
  beta.true_new <- beta.true[c("interc", col_all_good), , drop=FALSE ]
  #print("new beta true")
  #print(beta.true_new)
  rownames(beta.true_new)<-c("interc", col_all_good)
  #beta.true
  #print(beta.true)
  # Response vector (2way)
  sigma.y <- 5
  y.3w <- data.frame(row.names=rownames(x.3w))
  set.seed(123)
  y.3w$obs <- beta.true$coeffs[1] + as.matrix(x.3w)%*%as.vector(beta.true$coeffs)[-1] + rnorm(nrow(y.3w), 0, sigma.y)
  y.3w$true <- beta.true$coeffs[1] + as.matrix(x.3w)%*%as.vector(beta.true$coeffs)[-1]
  #print(colnames(x.3w))
  x.3w_new<-x.3w[,col_all_good ,drop=FALSE]
  
  if(all(rownames(beta.true_new)[-1]==colnames(x.3w_new)) ==TRUE)
  {print("Data loaded properly")}
  print("SNR:")
  print(var(as.vector(y.3w$true))/ (sigma.y^2) )
  #print("newcolnames")
  #print(colnames(x.3w_new))
  return(list('X'=as.matrix(x.3w_new), 'y'=y.3w, 'beta'=beta.true_new))
}




create_hier_dataset_paper_many_main_vary_interactions_2way<-function(magnitude_scale=3){
  set.seed(123)
  x.3w <- dummy.matrix(NF=3, NL=c(7,6,5))
  l1=6
  l2=5
  l3=4
  col_mains<-colnames(x.3w)[c(1:(l1+l2+l3))]
  #print(col_mains)
  col_psi<-colnames(x.3w)[c ( (l1+l2+l3+ l1*l2+ l2*l3 + l1*l3 + 1): (l1+l2+l3+ l1*l2+ l2*l3+ +l1*l3+ l1*l2*l3) )]
  col_theta_good<-c()
  for (i in c(1:l1)) {
    for (j in c(1:l2)) {
      # Create the string and append to the vector
      col_theta_good <- c(col_theta_good, paste0("A.", i, ":B.", j))
    }
    for(k in c(1:l3))
    {col_theta_good <- c(col_theta_good, paste0("A.", i, ":C.", k))}
  }
  for (j in c(1:l2))
  {for(k in c(1:l3))
  {col_theta_good <- c(col_theta_good, paste0("B.", j, ":C.", k))}
  }
  #print(col_theta_good)
  col_all_good=c(col_mains,col_theta_good,col_psi)
  #print(col_all_good)
  # Hierarchical Coefficients (2way)
  p.3w <- ncol(x.3w)
  n.3w <- p.3w + 1
  beta.min <- 1
  beta.max <- 10
  beta.true <- data.frame(rep(0, n.3w))
  rownames(beta.true) <- c("interc", colnames(x.3w))
  colnames(beta.true) <- c("coeffs")
  set.seed(123)
  beta.true$coeffs <- runif(n.3w, beta.min, beta.max)*sample(c(1,-1),size=n.3w,replace=TRUE)
  
  
  levs.true <- c("interc","A.1", "A.2", "A.3","A.4",  "B.1", "B.2","B.3", "C.1", "C.2","C.3",
                 "A.1:B.1","A.1:B.2","A.2:B.1","A.2:B.2",
                 "A.1:C.1","A.1:C.2", "A.2:C.1","A.2:C.2",
                 "B.1:C.1","B.3:C.3")
  levs.high<-c("A.1:B.1","A.1:B.2","A.2:B.1","A.2:B.2",
               "A.1:C.1","A.1:C.2", "A.2:C.1","A.2:C.2",
               "B.1:C.1","B.3:C.3")
  
  beta.true$coeffs[which(!is.element(rownames(beta.true),levs.true))] <- 0
  beta.true$coeffs[which(is.element(rownames(beta.true),levs.high))]<- magnitude_scale*beta.true$coeffs[which(is.element(rownames(beta.true),levs.high))]
  beta.true_new <- beta.true[c("interc", col_all_good), , drop=FALSE ]
  #print("new beta true")
  #print(beta.true_new)
  rownames(beta.true_new)<-c("interc", col_all_good)
  #beta.true
  #print(beta.true)
  # Response vector (2way)
  sigma.y <- 5
  y.3w <- data.frame(row.names=rownames(x.3w))
  set.seed(123)
  y.3w$obs <- beta.true$coeffs[1] + as.matrix(x.3w)%*%as.vector(beta.true$coeffs)[-1] + rnorm(nrow(y.3w), 0, sigma.y)
  y.3w$true <- beta.true$coeffs[1] + as.matrix(x.3w)%*%as.vector(beta.true$coeffs)[-1]
  #print(colnames(x.3w))
  x.3w_new<-x.3w[,col_all_good ,drop=FALSE]
  
  if(all(rownames(beta.true_new)[-1]==colnames(x.3w_new)) ==TRUE)
  {print("Data loaded properly")}
  print("SNR:")
  print(var(as.vector(y.3w$true))/ (sigma.y^2) )
  #print("newcolnames")
  #print(colnames(x.3w_new))
  return(list('X'=as.matrix(x.3w_new), 'y'=y.3w, 'beta'=beta.true_new))
}



set.seed(455)
create_basic_dataset<-function(){
  set.seed(455)
  x.3w <- dummy.matrix(NF=3, NL=c(9,9,5))
  l1=8
  l2=8
  l3=4
  col_mains<-colnames(x.3w)[c(1:(l1+l2+l3))]
  #print(col_mains)
  col_psi<-colnames(x.3w)[c ( (l1+l2+l3+ l1*l2+ l2*l3 + l1*l3 + 1): (l1+l2+l3+ l1*l2+ l2*l3+ +l1*l3+ l1*l2*l3) )]
  col_theta_good<-c()
  for (i in c(1:l1)) {
    for (j in c(1:l2)) {
      # Create the string and append to the vector
      col_theta_good <- c(col_theta_good, paste0("A.", i, ":B.", j))
    }
    for(k in c(1:l3))
    {col_theta_good <- c(col_theta_good, paste0("A.", i, ":C.", k))}
  }
  for (j in c(1:l2))
  {for(k in c(1:l3))
  {col_theta_good <- c(col_theta_good, paste0("B.", j, ":C.", k))}
  }
  #print(col_theta_good)
  col_all_good=c(col_mains,col_theta_good,col_psi)
  #print(col_all_good)
  # Hierarchical Coefficients (2way)
  p.3w <- ncol(x.3w)
  n.3w <- p.3w + 1
  beta.min <- 1
  beta.max <- 10
  beta.true <- data.frame(rep(0, n.3w))
  rownames(beta.true) <- c("interc", colnames(x.3w))
  colnames(beta.true) <- c("coeffs")
  set.seed(123)
  beta.true$coeffs <- runif(n.3w, beta.min, beta.max)*sample(c(1,-1),size=n.3w,replace=TRUE)
  
  
  levs.true <- c("A.1", "A.2", "A.3","A.4",  "B.1", "B.2","B.3", "C.1", "C.2","C.3",
                 "A.1:B.1","A.1:B.2","A.2:B.1","A.2:B.2",
                 "A.1:C.1","A.1:C.2", "A.2:C.1","A.2:C.2",
                 "B.1:C.1","B.1:C.2", "B.2:C.1", 
                 "A.1:B.1:C.1","A.1:B.1:C.2","A.2:B.1:C.1", "A.2:B.2:C.1")
  
  beta.true$coeffs[which(!is.element(rownames(beta.true),levs.true))] <- 0
  beta.true_new <- beta.true[c("interc", col_all_good), , drop=FALSE ]
  #print("new beta true")
  #print(beta.true_new)
  rownames(beta.true_new)<-c("interc", col_all_good)
  #beta.true
  #print(beta.true)
  # Response vector (2way)
  sigma.y <- 5
  y.3w <- data.frame(row.names=rownames(x.3w))
  set.seed(123)
  y.3w$obs <- beta.true$coeffs[1] + as.matrix(x.3w)%*%as.vector(beta.true$coeffs)[-1] + rnorm(nrow(y.3w), 0, sigma.y)
  y.3w$true <- beta.true$coeffs[1] + as.matrix(x.3w)%*%as.vector(beta.true$coeffs)[-1]
  #print(colnames(x.3w))
  x.3w_new<-x.3w[,col_all_good ,drop=FALSE]
  
  if(all(rownames(beta.true_new)[-1]==colnames(x.3w_new)) ==TRUE)
  {print("Data loaded properly")}
  print("SNR:")
  print(var(as.vector(y.3w$true))/ (sigma.y^2) )
  #print("newcolnames")
  #print(colnames(x.3w_new))
  return(list('X'=as.matrix(x.3w_new), 'y'=y.3w, 'beta'=beta.true_new))
}







