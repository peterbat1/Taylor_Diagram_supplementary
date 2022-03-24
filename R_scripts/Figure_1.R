# Taylor Diagrams for ENMs: Generate a simple trivariate Normal data set and
# plot a Taylor Diagram. This diagram appears as Figure 1 in the manuscript.
#
# The method used to generate a multivariate Mormal data set drew inspiration from the following source:
# https://blog.revolutionanalytics.com/2016/08/simulating-form-the-bivariate-normal-distribution-in-r-1.html
# last accessed 2022-03-10
#
#
# Peter D. Wilson
# Adjunct Fellow
# School of Natural Sciences
# Faculty of Science and Engineering
# Macquarie University, North Ryde, NSW, Australia 2109
#
# 2022-03-22 

library(MASS)
library(ggplot2)

##### Please change the path to the Taylor Diagram script to suit your system:
source("/home/peterw/Data_and_Projects/Personal Projects/Taylor Diagrams and ENMs/R-scripts/taylor_diagram_ggplot2_ver2.R")

N <- 200 # Number of random samples
set.seed(123)

# Target parameters for component univariate normal distributions
rho12 <- 0.6
rho13 <- -0.7
rho23 <- -0.2
mu1 <- 1; s1 <- 2
mu2 <- 1; s2 <- 4
mu3 <- 1; s3 <- 4

# Parameters for TRIVARIATE normal distribution
mu <- c(mu1, mu2, mu3) # Mean
sigma <- matrix(c(s1^2, s1*s2*rho12, s1*s3*rho13,
                   s1*s2*rho12, s2^2, s2*s3*rho23,
                   s1*s3*rho13, s2*s3*rho23, s3^2),
                 3, 3, byrow = TRUE) # Covariance matrix

bvn_data <- mvrnorm(N, mu = mu, Sigma = sigma ) # from MASS package
colnames(bvn_data) <- c("ref", "m1", "m2")

# OBSERVED summary statistics as a quality/sanity check:
mean_1 <- mean(bvn_data[, 1])
mean_2 <- mean(bvn_data[, 2])
mean_3 <- mean(bvn_data[, 3])

rmsd1_2 <- sqrt(sum(((bvn_data[, 2] - mean_2) - (bvn_data[, 1] - mean_1))^2)/N)
rmsd1_3 <- sqrt(sum(((bvn_data[, 3] - mean_2) - (bvn_data[, 1] - mean_1))^2)/N)

correl1_2 <- cor(bvn_data[, 1], bvn_data[, 2])
correl1_3 <- cor(bvn_data[, 1], bvn_data[, 3])

sd_1 <- sd(bvn_data[, 1])
sd_2 <- sd(bvn_data[, 2])
sd_3 <- sd(bvn_data[, 3])

plot(bvn_data[, c(1, 2)], main = "Positive correlation contrast")
plot(bvn_data[, c(1, 3)], main = "Negative correlation contrast")

# Plot the diagram
thisPlot <- taylor_diagram(data = bvn_data,
                           plot_type = "full",
                           show_labels = TRUE,
                           model_labels = c("Ref", "Model 1", "Model 2"))

# Please change the path/filename to suit your system:
ggsave("/home/peterw/Data_and_Projects/Personal Projects/Taylor Diagrams and ENMs/Results/Figure_1_new.png",
       thisPlot,
       device = "png",
       units = "mm",
       width = 160,
       height =160,
       bg = "white")

