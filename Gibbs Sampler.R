
library(coda)
y1<-read.table("http://www2.stat.duke.edu/~pdh10/FCBS/Exercises/school1.dat")
y1$school<-1
y2<-read.table("http://www2.stat.duke.edu/~pdh10/FCBS/Exercises/school2.dat")
y2$school<-2
y3<-read.table("http://www2.stat.duke.edu/~pdh10/FCBS/Exercises/school3.dat")
y3$school<-3
y4<-read.table("http://www2.stat.duke.edu/~pdh10/FCBS/Exercises/school4.dat")
y4$school<-4
y5<-read.table("http://www2.stat.duke.edu/~pdh10/FCBS/Exercises/school5.dat")
y5$school<-5
y6<-read.table("http://www2.stat.duke.edu/~pdh10/FCBS/Exercises/school6.dat")
y6$school<-6
y7<-read.table("http://www2.stat.duke.edu/~pdh10/FCBS/Exercises/school7.dat")
y7$school<-7
y8<-read.table("http://www2.stat.duke.edu/~pdh10/FCBS/Exercises/school8.dat")
y8$school<-8
y<-rbind(y1,y2,y3,y4,y5,y6,y7,y8)
colnames(y)[1] <-"hour"
Y <- y[c("school","hour")]

# prior parameters
nu0 <- 2; sigma20 <- 15
eta0 <- 2; tau20 <- 10
mu0 <- 7; gamma20 <- 5
# starting values
J <- length(unique(Y[ ,1]))
n.vec <- var.vec <- ybars <- rep(NA,J)
for(j in 1:J){
  ybars[j] <- mean(Y[Y[ ,1]==j, 2])
  var.vec[j] <- var(Y[Y[ ,1]==j, 2])
  n.vec[j] <- sum(Y[ ,1]==j)
}
muj.vec <- ybars; sigma2 <- mean(var.vec)
mu <- mean(muj.vec); tau2 <- var(muj.vec)
# Gibbs sampler
set.seed(108)
S <- 10000; burnin <- 0.3*S
muj.mat <- matrix(nrow=S, ncol=J)
other.mat <- matrix(nrow=S, ncol=3)
for(s in 1:(S + burnin)){
  # sample new values of the mujs
  for (j in 1:J){
    muj.var <- 1/(n.vec[j]/sigma2 + 1/tau2)
    muj.mean <- muj.var*(ybars[j]*n.vec [j]/sigma2 + mu/tau2)
    muj.vec[j] <- rnorm(1, muj.mean, sqrt(muj.var))
  }
  #sample new value of sigma2
  mu_n <- nu0 + sum(n.vec)
  gamma2_n <- nu0*sigma20
  for (j in 1:J){
    gamma2_n <- gamma2_n + sum((Y[Y[ ,1]==j ,2] - muj.vec[j])^2)
  }
  sigma2 <- 1/rgamma(1, mu_n/2, gamma2_n/2)
  #sample a new value of mu
  mu.var <- 1/(J/tau2 + 1/gamma20)
  mu.mean <- mu.var*(J*mean(muj.vec)/tau2 + mu0/gamma20)
  mu <- rnorm(1, mu.mean, sqrt(mu.var))
  # sample a new value o f tau2
  eta_n <- eta0 + J
  gamma2_n <- eta0*tau20 + sum((muj.vec-mu)^2)
  tau2 <- 1/rgamma(1, eta_n/2, gamma2_n/2)
  # store results
  if(s > burnin){
    muj.mat[s-burnin,] <- muj.vec
    other.mat[s-burnin,] <- c(sigma2, mu, tau2)
  }
}

mu.sample <- other.mat[,2]
winsize <- 500
mu.forplot <- matrix(NA, nrow=winsize, ncol=S/winsize)
for(w in 1:(S/winsize)){
  id1 <- (w-1)*winsize + 1
  id2 <- w*winsize
  mu.forplot[,w] = mu.sample[id1:id2]
}
par(las=2)
boxplot(mu.forplot,xlab="iteration",ylab="mu",main="Stationary Plot")

#e
# Posterior mean for sigma2
sigma2.sample <- other.mat[,1]
meansig <- mean(sigma2.sample)
cat("Posterior mean for sigma2", meansig)

# 95% quantile-based interval
qhighersig <- quantile(sigma2.sample,0.975)
qlowersig <- quantile(sigma2.sample,0.025)
# Posterior mean for tau2
tau2.sample <- other.mat[,3]
meantau<- mean(tau2.sample)
cat("Posterior mean for tau2", meantau)

# 95% quantile-based interval
qhighertau <- quantile(tau2.sample,0.975)
qlowertau <- quantile(tau2.sample,0.025)


