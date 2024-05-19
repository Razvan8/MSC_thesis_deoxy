# Clear all
rm(list=ls())
gc()
try(dev.off(dev.list()["RStudioGD"]),silent=TRUE)
try(dev.off(),silent=TRUE)
cat("\014")

# Memory
memory.limit(64000)



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
beta[128]=-2  #A1B9C3


noise<-rnorm(dim(x.3w)[1],0, 0.03 )
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
  print(noise_variance)
  # Generate the response
  Y <- signal + epsilon
  return(Y)
}

# Generate data with different SNR values
Y_high_snr <- generate_y_with_snr(X, beta, snr = 10)  # High SNR
Y_moderate_snr <- generate_y_with_snr(X, beta, snr = 5)  # Moderate SNR
Y_low_snr <- generate_y_with_snr(X, beta, snr = 0.5)  # Low SNR

# Check the variances and SNRs
signal_variance <- var(as.vector(X %*% beta))
high_snr <- signal_variance / var(Y_high_snr - X %*% beta)
moderate_snr <- signal_variance / var(Y_moderate_snr - X %*% beta)
low_snr <- signal_variance / var(Y_low_snr - X %*% beta)

print(paste("High SNR:", high_snr))
print(paste("Moderate SNR:", moderate_snr))
print(paste("Low SNR:", low_snr))




return(list('X'=X, 'y'=y))

}

#data=create_sparse_hier_dataset()
#X<-data$X
#y<-data$y
#print(dim(X))
#print(dim(y))




##Hierarchical datasets big mains smaller 2way smaller 3way




##Hierarchical dataset small mains, beigger 2ay 3way



##Hierarchical dataset no 3 way ()


##Non-Hierarchical dataset

