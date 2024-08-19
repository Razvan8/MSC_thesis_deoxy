library(glmnet)


set.seed(123)

col1 <- c(rep(1, 40), rep(0, 40), rep(-1, 40))
col2 <- c(rep(0, 40), rep(1, 40), rep(-1, 40))
mat <- cbind(col1, col2)

eps1 <- rnorm(n = nrow(mat), mean = 0, sd = 2)
y1 <- 3 * mat[,1] - 3 * mat[,2] + eps1


lasso_model1 <- glmnet(mat, y1, alpha = 1, lambda = 0.5)


coefficients1 <- coef(lasso_model1)


print("Coefficients for Case 1 (y = 3 * col1 - 3 * col2 + eps(0,1)):")
print(coefficients1)





eps2 <- rnorm(n = nrow(mat), mean = 0, sd = 1)
y2 <- 3 * mat[,1] + 0 * mat[,2] + eps2


lasso_model2 <- glmnet(mat, y2, alpha = 1, lambda = 0.5)


coefficients2 <- coef(lasso_model2)






# Define the function


print("Coefficients for Case 2 (y = 3 * col1 + 0 * col2 + eps(0,1)):")
print(coefficients2)







# Define the function
g.link <- function(x) {
  tt <- apply(as.matrix(x), 1, FUN = function(v) min(max(v, 0.001), 0.999))
  g <- 3.5 * tan(pi * (2 * tt - 1) / 2)
  return(g)
}

# Generate the data
x <- seq(0.06, 0.94, length.out = 100)
y <- g.link(x)

# Define the function
g.link <- function(x) {
  tt <- apply(as.matrix(x), 1, FUN = function(v) min(max(v, 0.001), 0.999))
  g <- 3.5 * tan(pi * (2 * tt - 1) / 2)
  return(g)
}

# Load the necessary library
library(ggplot2)

# Define the function
g.link <- function(x) {
  tt <- apply(as.matrix(x), 1, FUN = function(v) min(max(v, 0.001), 0.999))
  g <- 3.5 * tan(pi * (2 * tt - 1) / 2)
  return(g)
}

# Generate the data
x <- seq(0.07, 0.93, length.out = 1000)
y <- g.link(x)
data <- data.frame(x, y)



element_textbox <- function(...) {
  el <- element_text(...)
  class(el) <- c("element_textbox", class(el))
  el
}


# Create the plot using ggplot2
ggplot(data, aes(x, y)) +
  geom_line(color = "blue") +
  ggtitle("Link function") +
  xlab("Expected Yield") +
  ylab("Systematic Component")+
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0, 1)) +
  theme(aspect.ratio=1)+
  theme(text = element_text(size = 12)) +
  theme(plot.title = element_textbox(hjust = 0.5, margin = margin(t = 5, b = 5))) 


