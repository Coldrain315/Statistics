library(coda)
data <- read.table("http://www2.stat.duke.edu/~pdh10/FCBS/Exercises/msparrownest.dat")
# print(data)
Y <- data$V1
wingspan <- data$V2
X <- model.matrix(~wingspan)
n <- nrow(X)
p <- ncol(X)
df <- data.frame(X)
df['Y'] <- Y
model_glm <- glm(Y~wingspan, data = df, family = "binomial")
summary(model_glm)
alpha_glm <- model_glm$coefficients[1]
beta_glm <- model_glm$coefficients[2]
predictions <- predict(model_glm, newdata=df, type="response")
sigma2_glm <- var(predictions)
sigma_0 <- n*sigma2_glm*solve(t(X)%*%X)
print(sigma_0)
mu_0 <- matrix(model_glm$coefficients)

# Metropolis
library(Rlab)
library(coda)
library(rsq)
library(mvtnorm)
#Set paramters for proposal density
c <- 0.5
delta <- 10 #use it to tune acceptance ratio
var_prop <- delta*var(log(Y+c))*solve(t(X)%*%X)
#Initial values for sampler
beta <- mu_0
#First set number of iterations and burn-in, then set seed
n_iter <- 10000
burn_in <- 0.3*n_iter
thin <- 5
set.seed (123)
#Set counter for acceptances
accept_counter <- 0
#Set null matrices to save samples
BETA <- matrix(0,nrow=n_iter,ncol=p)

#Now, to the sampler
for(s in 1:(n_iter+burn_in)){
  #generate proposal
  beta_star <- t(rmvnorm(1,beta,var_prop))
  #compute acceptance ratio/probability
  #do so on log scale because r can be numerically unstable
  theta <- exp(X%*%beta_star)/(1+exp(X%*%beta_star))
  log_r <- sum(dbern(Y,prob=theta),log=T) +
    dmvnorm(c(beta_star),mu_0,sigma_0,log=T) -
    sum(dbern(Y,theta),log=T) -
    dmvnorm(c(beta),mu_0,sigma_0,log=T)
  if(log(runif(1)) < log_r){
    accept_counter <- accept_counter + 1
    beta <- beta_star
  }
  if(s > burn_in){
    BETA[(s-burn_in),] <- beta 
  }
}
#Check acceptance rate
accept_counter/(n_iter+burn_in)

#thinning
sample_thin <- seq(1,n_iter,by=thin)
BETA_thinned <- BETA[sample_thin,]
library(coda)
alpha_chain <- BETA_thinned[,1]
beta_chain <- BETA_thinned[,2]
effectiveSize(mcmc(alpha_chain))
effectiveSize(mcmc(beta_chain))
plot(mcmc(BETA_thinned))