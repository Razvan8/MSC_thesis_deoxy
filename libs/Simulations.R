####Sparse hierarchical simulations####

libs_path<-file.path("..","libs")
source(file.path(libs_path,'Create_synthetic_datasets.R'))


data=create_sparse_hier_dataset()
X<-data$X
y<-data$y
print(dim(X))
print(dim(y))