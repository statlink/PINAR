
## test TNAR model bootstrap

library(igraph)
library(copula)
library(nloptr)

T1 <- 1000 # time
NN <- 20   # number of nodes

W2 <- adja(NN, 2, 0.3)   # network

# Lag 1, no covariates

y2 <- poisson.MODpq.tar(c(0.2,0.25, 0.20), W2, gamma=1, a=c(0.2, 0.2, 0.2), p=1, d=1, Z=NULL,
                         Time=T1, N=NN, copula="gaussian", rho=0.5, df=1)$y # data

x02 <- ols.narpq(N=NN, TT=T1, y2, W2, p=1, Z=NULL) 
x02      ## ols as starting values optimization

lin2 <- lin_estimnarpq(x02, NN, T1, y2, W2, p=1, Z=NULL)
lin2    ## estimates of linear model

#Find interval for Brent algorithm computed as mean(20%-quantiles) - mean(80%-quantiles)
# in the distribution of W%*%y[,t]

Z2=W2%*%y2
qq <- matrix(0, 2, nrow(Z2))
for (k in 1:nrow(Z2)){
  qq[,k] <- as.matrix(quantile(Z2[k,], prob=c(0.20, 0.80)))
}

mean(qq[1,])
mean(qq[2,])

g_lm <- global_optimise_LM(f=LM_gamma_tpnarpq, L=mean(qq[1,]), U=mean(qq[2,]), I=2, lin2$coeflin, NN, T1,
                           y2, W2, p=1, d=1, Z=NULL)
g_lm  ## optimal value gamma, sup test and intervals



boot2 <- score_test_tarpq_j(g_lm$supLM, lin2$coeflin, NN, T1, y2, W2, p=1, d=1,
                             Z=NULL, J=199, L=mean(qq[1,]), U=mean(qq[2,]))
boot2      ## test with bootstrap, by optimizing the test at each j 
boot2$cpJ


# Lag 2 and covariates

Z <- matrix(rexp(NN*2),NN,2) # positive covariates

y3 <- poisson.MODpq.tar(c(1,0.2, 0.1, 0.1, 0.1,0.2,0.3), W2, gamma=1, a=c(0.1, 0.15, 0.1, 0.1, 0.1), p=2, d=1, Z=Z,
                         Time=T1, N=NN, copula="gaussian", rho=0.5, df=1)$y # data

x03 <- ols.narpq(N=NN, TT=T1, y3, W2, p=2, Z=Z) 
x03      ## ols as starting values optimization

lin3 <- lin_estimnarpq(x03, NN, T1, y3, W2, p=2, Z=Z)
lin3    ## estimates of linear model

# find quantile interval

Z3=W2%*%y3
qq3 <- matrix(0, 2, nrow(Z3))
for (k in 1:nrow(Z3)){
  qq3[,k] <- as.matrix(quantile(Z3[k,], prob=c(0.20, 0.80)))
}

mean(qq3[1,])
mean(qq3[2,])

g_lm3 <- global_optimise_LM(f=LM_gamma_tpnarpq, L=mean(qq3[1,]), U=mean(qq3[2,]), I=2, lin3$coeflin, NN, T1,
                            y3, W2, p=2, d=1, Z=Z)
g_lm3  ## optimal value gamma, sup test and intervals


# d=2

g_lm3.2 <- global_optimise_LM(f=LM_gamma_tpnarpq, L=mean(qq3[1,]), U=mean(qq3[2,]),
                              I=2, lin3$coeflin, NN, T1, y3, W2, p=2, d=2, Z=Z)
g_lm3.2



boot3.3 <- score_test_tarpq_j(g_lm3$supLM, lin3$coeflin, NN, T1, y3, W2, p=2, d=1,
                            Z=Z, J=199, L=mean(qq3[1,]), U=mean(qq3[2,]))
boot3.3      ## test with bootstrap, by optimizing the test at each j 
boot3.3$pJ


boot3.4 <- score_test_tarpq_j(g_lm3.2$supLM, lin3$coeflin, NN, T1, y3, W2, p=2, d=2,
                              Z=Z, J=199, L=mean(qq3[1,]), U=mean(qq3[2,]))
boot3.4      ## test with bootstrap, by optimizing the test at each j 
boot3.4$pJ

######## under H0 #################

# Lag 1, no covariates

y0 <- poisson.MODpq(c(1,0.2, 0.2), W2, p=1, Z=NULL,
                    Time=T1, N=NN, copula="gaussian", rho=0.5, df=1)$y # data

x00 <- ols.narpq(N=NN, TT=T1, y0, W2, p=1, Z=NULL) 
x00      ## ols as starting values optimization

lin0 <- lin_estimnarpq(x00, NN, T1, y0, W2, p=1, Z=NULL)
lin0    ## estimates of linear model

# find quantiles 

Z0=W2%*%y0
qq0 <- matrix(0, 2, nrow(Z0))
for (k in 1:nrow(Z0)){
  qq0[,k] <- as.matrix(quantile(Z0[k,], prob=c(0.20, 0.80)))
}

mean(qq0[1,])
mean(qq0[2,])

g_lm0 <- global_optimise_LM(f=LM_gamma_pnarpq, L=mean(qq0[1,]), U=mean(qq0[2,]),
                            I=2, lin0$coeflin, NN, T1, y0, W2, p=1, d=1, Z=NULL)
g_lm0  ## optimal value gamma, sup test and intervals


boot0 <- score_test_tarpq_j(g_lm0$supLM, lin0$coeflin, NN, T1, y0, W2, p=1, d=1,
                             Z=NULL, J=199, L=mean(qq0[1,]), U=mean(qq0[2,]))
boot0      ## test with bootstrap, by optimizing the test at each j 
boot0$cpJ

# Lag 2 and covariates

y4 <- poisson.MODpq(c(1,0.2, 0.2, 0.1, 0.1,0.2,0.3), W2, p=2, Z=Z,
                    Time=T1, N=NN, copula="gaussian", rho=0.5, df=1)$y # data

x04 <- ols.narpq(N=NN, TT=T1, y4, W2, p=2, Z=Z) 
x04      ## ols as starting values optimization

lin4 <- lin_estimnarpq(x04, NN, T1, y4, W2, p=2, Z=Z)
lin4    ## estimates of linear model

# find quantiles 

Z4=W2%*%y4
qq4 <- matrix(0, 2, nrow(Z4))
for (k in 1:nrow(Z4)){
  qq4[,k] <- as.matrix(quantile(Z4[k,], prob=c(0.20, 0.80)))
}

mean(qq4[1,])
mean(qq4[2,])


g_lm4 <- global_optimise_LM(f=LM_gamma_tpnarpq, L=mean(qq4[1,]), U=mean(qq4[2,]),
                            I=2, lin4$coeflin, NN, T1, y4, W2, p=2, d=1, Z=Z)
g_lm4  ## optimal value gamma, sup test and intervals


# d=2

g_lm4.2 <- global_optimise_LM(f=LM_gamma_tpnarpq, L=mean(qq4[1,]), U=mean(qq4[2,]),
                              I=2, lin4$coeflin, NN, T1, y4, W2, p=2, d=2, Z=Z)
g_lm4.2


boot5 <- score_test_tarpq_j(g_lm4$supLM, lin4$coeflin, NN, T1, y4, W2, p=2, d=1,
                            Z=Z, J=199, L=mean(qq4[1,]), U=mean(qq4[2,]))
boot5      ## test with bootstrap, by optimizing the test at each j 
boot5$pJ

boot5.2 <- score_test_tarpq_j(g_lm4.2$supLM, lin4$coeflin, NN, T1, y4, W2, p=2, d=2,
                              Z=Z, J=199, L=mean(qq4[1,]), U=mean(qq4[2,]))
boot5.2      ## test with bootstrap, by optimizing the test at each j 
boot5.2$pJ







## ignore all these
#
#
# u0 <- seq(from=mean(qq0[1,]), to=mean(qq0[2,]), length.out = 30)
# g_u0 <- vector()
# 
# for (k in 1:30) {
#   g_u0[k] <- LM_gamma_tpnarpq(u0[k], lin0$coeflin, NN, T1, y0, W2, p=1, d=1, Z=NULL)
# }
# 
# plot(g_u0 ~ u0) 
# g_u0        
#
/
# u3.2 <- seq(from=mean(qq3[1,]), to=mean(qq3[2,]), length.out = 30)
# g_u3.2 <- vector()
# 
# for (k in 1:30) {
#   g_u3.2[k] <- LM_gamma_tpnarpq(u3.2[k], lin3$coeflin, NN, T1, y3, W2, p=2, d=2, Z=Z)
# }
# 
# plot(g_u3.2 ~ u3.2) 
# g_u3.2          
#
# 
# boot3 <- score_test_tarpq_j(g_lm3$supLM, lin3$coeflin, NN, T1, y3, W2, p=2, d=1,
#                             Z=Z, J=199, L=mean(qq3[1,]), U=mean(qq3[2,]))
# boot3      ## test with bootstrap, by optimizing the test at each j 
# boot3$pJ
# 
# 
# boot3.2 <- score_test_tarpq_j(g_lm3.2$supLM, lin3$coeflin, NN, T1, y3, W2, p=2, d=2,
#                               Z=Z, J=199, L=mean(qq3[1,]), U=mean(qq3[2,]))
# boot3.2      ## test with bootstrap, by optimizing the test at each j 
# boot3.2$pJ
#
#
# u <- seq(from=0.5, to=mean(qq[2,]), length.out = 20)
# g_u <- vector()
# 
# for (k in 1:20) {
#   g_u[k] <- LM_gamma_tpnarpq(u[k], lin2$coeflin, NN, T1, y2, W2, p=1, d=1, Z=NULL)
# }
# 
# plot(g_u ~ u) # The test is a step function. Although there is not a unique value
# g_u           # of gamma for the maximum of the test the optimizer seems to work fine
#               # if it does not just change the value of L and U in the optimizer
# quantile(Z2, prob=c(0.1, 0.2, 0.25, 0.5, 0.75, 0.8, 0.9))
# 
# boot2q <- score_test_tarpq_j(g_lm$supLM, lin2$coeflin, NN, T1, y2, W2, p=1, d=1,
#                             Z=NULL, J=199, L=0.3, U=1.5)
# boot2q     ## test with bootstrap, by optimizing the test at each j
# boot2q$pJ
# u3 <- seq(from=mean(qq3[1,]), to=mean(qq3[2,]), length.out = 30)
# g_u3 <- vector()
# 
# for (k in 1:30) {
#   g_u3[k] <- LM_gamma_tpnarpq(u3[k], lin3$coeflin, NN, T1, y3, W2, p=2, d=1, Z=Z)
# }
# 
# plot(g_u3 ~ u3) 
# g_u3          