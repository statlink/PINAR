
## test STNAR model bootstrap

library(igraph)
library(copula)
library(nloptr)

T1 <- 1000 # time
NN <- 20   # number of nodes

W2 <- adja(NN, 2, 0.3)   # network

# Lag 1, no covariates

y2 <- poisson.MODpq.star(c(1,0.2, 0.2), W2, gamma=0.1, a=0.4, p=1, d=1, Z=NULL,
                   Time=T1, N=NN, copula="gaussian", rho=0.5, df=1)$y # data

x02 <- ols.narpq(N=NN, TT=T1, y2, W2, p=1, Z=NULL) 
x02      ## ols as starting values optimization

lin2 <- lin_estimnarpq(x02, NN, T1, y2, W2, p=1, Z=NULL)
lin2    ## estimates of linear model

#Find interval for Brent algorithm computed such that 0.1 <= exp(-gamma* X^2) <=0.9
# where X is the overall mean of W%*%y

Z2=W2%*%y2
exp(-0.07*(mean(Z2))^2)
exp(-1.45*(mean(Z2))^2)

g_lm <- global_optimise_LM(f=LM_gamma_stpnarpq, L=0.07, U=1.45, I=6, lin2$coeflin, NN, T1,
                           y2, W2, p=1, d=1, Z=NULL)
g_lm  ## optimal value gamma, sup test and intervals 
      ## (this also gives a guess of possible interval of values for DV below)

DV2 <- score_test_starpq_DV(lin2$coeflin, NN, T1, y2, W2, p=1, d=1, Z=NULL,
                            gamma_L=0.07, gamma_U=1.45, l=100) 
DV2     ## test vs STNAR with Davies bound


boot2.3 <- score_test_starpq_j(g_lm$supLM, lin2$coeflin, NN, T1, y2, W2, p=1, d=1,
                             Z=NULL, J=199,  L=0.07, U=1.45)
boot2.3      ## test with bootstrap, by optimizing the test at each j
boot2.3$cpJ


# Lag 2 and covariates

Z <- matrix(rexp(NN*2),NN,2) # positive covariates

y3 <- poisson.MODpq.star(c(1,0.1, 0.1, 0.1, 0.1,0.2,0.3), W2, gamma=0.1, a=c(0.4,0.1), p=2, d=1, Z=Z,
                         Time=T1, N=NN, copula="gaussian", rho=0.5, df=1)$y # data

x03 <- ols.narpq(N=NN, TT=T1, y3, W2, p=2, Z=Z) 
x03      ## ols as starting values optimization

lin3 <- lin_estimnarpq(x03, NN, T1, y3, W2, p=2, Z=Z)
lin3    ## estimates of linear model

# find gamma interval

Z3=W2%*%y3
exp(-0.03*(mean(Z3))^2)
exp(-0.65*(mean(Z3))^2)

g_lm3 <- global_optimise_LM(f=LM_gamma_stpnarpq, L=0.03, U=0.65, I=2, lin3$coeflin, NN, T1,
                           y3, W2, p=2, d=1, Z=Z)
g_lm3  ## optimal value gamma, sup test and intervals

DV3 <- score_test_starpq_DV(b=lin3$coeflin, N=NN, TT=T1, y=y3, W=W2, p=2, d=1, Z=Z,
                            gamma_L=0.03, gamma_U=0.65, l=100) 
DV3     ## test vs STNAR with Davies bound


# optimization with d=2

g_lm3.2 <- global_optimise_LM(f=LM_gamma_stpnarpq, L=0.03, U=0.65, I=2, lin3$coeflin, NN, T1,
                            y3, W2, p=2, d=2, Z=Z)
g_lm3.2

DV3.2 <- score_test_starpq_DV(lin3$coeflin, NN, T1, y3, W2, p=2, d=2, Z=Z,
                              gamma_L=0.03, gamma_U=0.65, l=100)
DV3.2     ## test vs STNAR with Davies bound


boot3 <- score_test_starpq_j(g_lm3$supLM, lin3$coeflin, NN, T1, y3, W2, p=2, d=1,
                             Z=Z, J=199, L=0.03, U=0.65)
boot3      ## test with bootstrap, by optimizing the test at each j 
boot3$pJ

boot3.2 <- score_test_starpq_j(g_lm3.2$supLM, lin3$coeflin, NN, T1, y3, W2, p=2, d=2,
                               Z=Z, J=199, L=0.03, U=0.65)
boot3.2      ## test with bootstrap, by optimizing the test at each j 
boot3.2$pJ



######## under H0 #################

# Lag 1, no covariates

y0 <- poisson.MODpq(c(1,0.2, 0.2), W2, p=1, Z=NULL,
                         Time=T1, N=NN, copula="gaussian", rho=0.5, df=1)$y # data

x00 <- ols.narpq(N=NN, TT=T1, y0, W2, p=1, Z=NULL) 
x00      ## ols as starting values optimization

lin0 <- lin_estimnarpq(x00, NN, T1, y0, W2, p=1, Z=NULL)
lin0    ## estimates of linear model

# find gamma interval

Z0=W2%*%y0
exp(-0.095*(mean(Z0))^2)
exp(-2.3*(mean(Z0))^2)

g_lm0 <- global_optimise_LM(f=LM_gamma_stpnarpq, L=0.095, U=2.6, I=2, lin0$coeflin,
                            NN, T1, y0, W2, p=1, d=1, Z=NULL)
g_lm0  ## optimal value gamma, sup test and intervals

DV0 <- score_test_starpq_DV(lin0$coeflin, NN, T1, y0, W2, p=1, d=1, Z=NULL,
                            gamma_L=0.095, gamma_U=2.3, l=100) 
DV0     ## test vs STNAR with Davies bound


boot0 <- score_test_starpq_j(g_lm0$supLM, lin0$coeflin, NN, T1, y0, W2, p=1, d=1,
                             Z=NULL, J=199, L=0.095, U=2.3)
boot0      ## test with bootstrap, by optimizing the test at each j 
boot0$cpJ

# Lag 2 and covariates

y4 <- poisson.MODpq(c(1,0.2, 0.2, 0.1, 0.1,0.2,0.3), W2, p=2, Z=Z,
                         Time=T1, N=NN, copula="gaussian", rho=0.5, df=1)$y # data

x04 <- ols.narpq(N=NN, TT=T1, y4, W2, p=2, Z=Z) 
x04      ## ols as starting values optimization

lin4 <- lin_estimnarpq(x04, NN, T1, y4, W2, p=2, Z=Z)
lin4    ## estimates of linear model

# find gamma interval

Z4=W2%*%y4
exp(-0.025*(mean(Z4))^2)
exp(-0.6*(mean(Z4))^2)

g_lm4 <- global_optimise_LM(f=LM_gamma_stpnarpq, L=0.025, U=0.6, I=2, lin4$coeflin,
                            NN, T1, y4, W2, p=2, d=1, Z=Z)
g_lm4  ## optimal value gamma, sup test and intervals

DV4 <- score_test_starpq_DV(b=lin4$coeflin, N=NN, TT=T1, y=y4, W=W2, p=2, d=1, Z=Z,
                            gamma_L=0.025, gamma_U=0.6, l=100) 
DV4     ## test vs STNAR with Davies bound

# optimization with d=2

g_lm4.2 <- global_optimise_LM(f=LM_gamma_stpnarpq, L=0.025, U=0.6, I=2, lin4$coeflin, NN, T1,
                              y4, W2, p=2, d=2, Z=Z)
g_lm4.2

DV4.2 <- score_test_starpq_DV(lin4$coeflin, NN, T1, y4, W2, p=2, d=2, Z=Z,
                              gamma_L=0.025, gamma_U=0.6, l=100)
DV4.2     ## test vs STNAR with Davies bound


boot5 <- score_test_starpq_j(g_lm4$supLM, lin4$coeflin, NN, T1, y4, W2, p=2, d=1,
                             Z=Z, J=199, L=0.025, U=0.6)
boot5      ## test with bootstrap, by optimizing the test at each j 
boot5$pJ

boot5.2 <- score_test_starpq_j(g_lm4.2$supLM, lin4$coeflin, NN, T1, y4, W2, p=2, d=2,
                               Z=Z, J=199, L=0.025, U=0.6)
boot5.2      ## test with bootstrap, by optimizing the test at each j 
boot5.2$pJ






# ignore all these
#
# boot2 <- score_test_starpq_j(g_lm$supLM, lin2$coeflin, NN, T1, y2, W2, p=1, d=1,
#                              Z=NULL, J=199,  L=0.1, U=1.2)
# boot2      ## test with bootstrap, by optimizing the test at each j
# boot2$cpJ
#
#
# boot2.2 <- score_test_starpq_j(g_lm$supLM, lin2$coeflin, NN, T1, y2, W2, p=1, d=1,
#                              Z=NULL, J=199,  L=0.0001, U=1000)
# boot2.2      ## test with bootstrap, by optimizing the test at each j 
# boot2.2$cpJ
#
#
# boot22 <- score_test_starpq_j2(lin2$coeflin, 20, 100, y2, W2, p=1, d=1, Z=NULL, J=200, U=1) 
# boot22      ## test with bootstrap, by using the same optimized gamma by the sample
# boot22$pJ
# boot4 <- score_test_starpq_j2(lin3$coeflin, 20, 100, y3, W2, p=2, d=1, Z=Z, J=200, U=1) 
# boot4      ## test with bootstrap, by using the same optimized gamma by the sample
# boot4$pJ
# 
# boot4.2 <- score_test_starpq_j2(lin3$coeflin, 20, 100, y3, W2, p=2, d=2, Z=Z, J=200, U=1) 
# boot4.2      ## test with bootstrap, by using the same optimized gamma by the sample
# boot4.2$pJ
