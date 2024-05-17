libs_path<-file.path("..","libs")
source(file.path(libs_path,'helper_functions.R'))
data<-load_deoxi_flourination()

print(max(data$yield))
hist(data$prob)

which(data$yield_corrected==100)
which(data$yield>100)

data


plot_yield_histo(data)


levels(data$a)
levels(data$b)
levels(data$s)
37*4*5

dim(data)