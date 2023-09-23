
## test

library(igraph)
library(copula)
library(nloptr)

W2 <- adja(20, 2, 0.3)   # network

# Lag 1, covariates

y2 <- poisson.MODpq.star(c(1,0.2, 0.2), W2, gamma=0.1, a=0.4, p=1, d=1, Z=NULL,
                   Time=100, N=20, copula="gaussian", rho=0.5, df=1)$y # data

x02 <- ols.narpq(N=20, TT=100, y2, W2, p=1, Z=NULL) 
x02      ## ols as starting values optimization

lin2 <- lin_estimnarpq(x02, 20, 100, y2, W2, p=1, Z=NULL)
lin2    ## estimates of linear model

g_lm <- global_optimise_LM(f=LM_gamma_stpnarpq, L=0.001, U=10, I=4, lin2$coeflin, 20, 100,
                   y2, W2, p=1, d=1, Z=NULL)
g_lm  ## optimal value gamma, sup test and intervals

DV2 <- score_test_starpq_DV(lin2$coeflin, 20, 100, y2, W2, p=1, d=1, Z=NULL,
                            gamma_L=0.001, gamma_U=10, l=100) 
DV2     ## test vs STNAR with Davies bound

boot2 <- score_test_starpq_j(g_lm$supLM, lin2$coeflin, 20, 100, y2, W2, p=1, d=1,
                             Z=NULL, J=199,  L=0.001, U=100)
boot2      ## test with bootstrap, by optimizing the test at each j 
boot2$cpJ


# Lag 2 and covariates

Z <- matrix(rexp(20*2),20,2)

y3 <- poisson.MODpq.star(c(1,0.1, 0.1, 0.1, 0.1,0.2,0.3), W2, gamma=0.1, a=c(0.4,0.1), p=2, d=1, Z=Z,
                         Time=100, N=20, copula="gaussian", rho=0.3, df=1)$y # data

x03 <- ols.narpq(N=20, TT=100, y3, W2, p=2, Z=Z) 
x03      ## ols as starting values optimization

lin3 <- lin_estimnarpq(x03, 20, 100, y3, W2, p=2, Z=Z)
lin3    ## estimates of linear model

g_lm3 <- global_optimise_LM(f=LM_gamma_stpnarpq, L=0.01, U=10, I=5, lin3$coeflin, 20, 100,
                           y3, W2, p=2, d=1, Z=Z)
g_lm3  ## optimal value gamma, sup test and intervals

DV3 <- score_test_starpq_DV(b=lin3$coeflin, N=20, TT=100, y=y3, W=W2, p=2, d=1, Z=Z,
                            gamma_L=0.01, gamma_U=10, l=100) 
DV3     ## test vs STNAR with Davies bound

# optimization with d=2

g_lm3.2 <- global_optimise_LM(f=LM_gamma_stpnarpq, L=0.01, U=10, I=4, lin3$coeflin, 20, 100,
                            y3, W2, p=2, d=2, Z=Z)
g_lm3.2

DV3.2 <- score_test_starpq_DV(lin3$coeflin, 20, 100, y3, W2, p=2, d=2, Z=Z,
                              gamma_L=0.01, gamma_U=10, l=100)
DV3.2     ## test vs STNAR with Davies bound
DV3.2$DV


boot3 <- score_test_starpq_j(g_lm3$supLM, lin3$coeflin, 20, 100, y3, W2, p=2, d=1, Z=Z, J=199, L=0.001, U=50)
boot3      ## test with bootstrap, by optimizing the test at each j 
boot3$pJ

boot3.2 <- score_test_starpq_j(g_lm3.2$supLM, lin3$coeflin, 20, 100, y3, W2, p=2, d=2, Z=Z, J=199, L=0.001, U=100)
boot3.2      ## test with bootstrap, by optimizing the test at each j 
boot3.2$pJ



######## under H0 #################

# Lag 1, no covariates

y0 <- poisson.MODpq(c(1,0.2, 0.2), W2, p=1, Z=NULL,
                         Time=100, N=20, copula="gaussian", rho=0.5, df=1)$y # data

x00 <- ols.narpq(N=20, TT=100, y0, W2, p=1, Z=NULL) 
x00      ## ols as starting values optimization

lin0 <- lin_estimnarpq(x00, 20, 100, y0, W2, p=1, Z=NULL)
lin0    ## estimates of linear model

g_lm0 <- global_optimise_LM(f=LM_gamma_stpnarpq, L=0.01, U=100, I=20, lin0$coeflin, 20, 100,
                           y0, W2, p=1, d=1, Z=NULL)
g_lm0  ## optimal value gamma, sup test and intervals

DV0 <- score_test_starpq_DV(lin0$coeflin, 20, 100, y0, W2, p=1, d=1, Z=NULL,
                            gamma_L=0.1, gamma_U=100, l=100) 
DV0     ## test vs STNAR with Davies bound


boot0 <- score_test_starpq_j(g_lm0$supLM, lin0$coeflin, 20, 100, y0, W2, p=1, d=1,
                             Z=NULL, J=199, L=0.001, U=100)
boot0      ## test with bootstrap, by optimizing the test at each j 
boot0$cpJ

# Lag 2 and covariates

y4 <- poisson.MODpq(c(1,0.1, 0.1, 0.1, 0.1,0.2,0.3), W2, p=2, Z=Z,
                         Time=100, N=20, copula="gaussian", rho=0.5, df=1)$y # data

x04 <- ols.narpq(N=20, TT=100, y4, W2, p=2, Z=Z) 
x04      ## ols as starting values optimization

lin4 <- lin_estimnarpq(x04, 20, 100, y4, W2, p=2, Z=Z)
lin4    ## estimates of linear model

g_lm4 <- global_optimise_LM(f=LM_gamma_stpnarpq, L=0.01, U=50, I=10, lin4$coeflin, 20, 100,
                            y4, W2, p=2, d=1, Z=Z)
g_lm4  ## optimal value gamma, sup test and intervals

DV4 <- score_test_starpq_DV(b=lin4$coeflin, N=20, TT=100, y=y4, W=W2, p=2, d=1, Z=Z,
                            gamma_L=0.01, gamma_U=50, l=100) 
DV4     ## test vs STNAR with Davies bound

# optimization with d=2

g_lm4.2 <- global_optimise_LM(f=LM_gamma_stpnarpq, L=0.01, U=100, I=10, lin4$coeflin, 20, 100,
                              y4, W2, p=2, d=2, Z=Z)
g_lm4.2

DV4.2 <- score_test_starpq_DV(lin4$coeflin, 20, 100, y4, W2, p=2, d=2, Z=Z,
                              gamma_L=0.01, gamma_U=20, l=100)
DV4.2     ## test vs STNAR with Davies bound
DV4.2$DV


boot5 <- score_test_starpq_j(g_lm4$supLM, lin4$coeflin, 20, 100, y4, W2, p=2, d=1, Z=Z, J=199, L=0.001, U=50)
boot5      ## test with bootstrap, by optimizing the test at each j 
boot5$pJ

boot5.2 <- score_test_starpq_j(g_lm4.2$supLM, lin4$coeflin, 20, 100, y4, W2, p=2, d=2, Z=Z, J=199, L=0.001, U=100)
boot5.2      ## test with bootstrap, by optimizing the test at each j 
boot5.2$pJ



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
