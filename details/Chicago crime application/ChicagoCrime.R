
################################################################################
## Chicago crime application for estimating linear and log-linear PNAR model
## and testing linearity of PNAR versus the following nonlinear models:
##    - Intercept Drift (ID)
##    - Smooth Transition (STNAR)
##    - Threshold (TNAR)
################################################################################

####Load required packages

library(GNAR)
library(igraph)
library(nloptr)
library(MASS)
library(ks)
library(car)
library(foreach)
library(readr)
library(Matrix)

set.seed(123)

# requires to load the functions in the following scripts:
#  - EstimatePnar1-2.R <-- functions to estimate linear PNAR model
#  - EstimateLogPnar1-2.R <-- functions to estimate log-linear PNAR model
#  - TestLinPNAR1vsID.R <-- functions to test linearity against ID nonlinear model
#  - TestLinPNAR1vsSTNAR.R <-- functions to test linearity against STNAR nonlinear model
#  - TestLinPNAR1vsTNAR.R <-- functions to test linearity against TNAR nonlinear model

# load data

crime <- read_csv("crime.csv")

net_crime <- readMM("neighborhood.mtx")

#generate adjacency matrix

gr_c <- graph_from_adjacency_matrix(net_crime, mode="undirected")
adj_c <- as_adjacency_matrix(gr_c, sparse = F) 
n_ic <- apply(adj_c, 1, sum)                  #out-degree
n_jc <- apply(adj_c, 2, sum)
diac <- diag(n_ic)
diac <- 1/diac
diac[which(diac==Inf)] <- 0
Wc <- diac%*%adj_c                     #row-normalized adjacency matrix
max(Wc)
min(Wc)

data <- as.matrix(crime[,-1])
W <- Wc


N <- nrow(W)      # number of nodes in the network
TT <- ncol(data)  # temporal sample size

## estimation of linear and log-linear PNAR models with 1 and 2 lags

# linear 

ols_chi <- ols.nar1(TT=TT, y=data, W=W)
ols_chi

estim_lin1_chi <- lin_estimnar1(x0=ols_chi, N=N, TT=TT, y=data, W=W) 
estim_lin1_chi

ols_2chi <- ols.nar2(N=N, TT=TT, y=data, W=W)
ols_2chi

estim_lin2_chi <- lin_estimnar2(x0=ols_2chi, N=N, TT=TT, y=data, W=W)
estim_lin2_chi

# log-linear 

ols_log_chi <- ols.nar1.log(N=N, TT=TT, y=data, W=W)
ols_log_chi

estim_log1_chi <- log_lin_estimnar1(x0=ols_log_chi, N=N, TT=TT, y=data, W=W) 
estim_log1_chi

ols_log_2chi <- ols.nar2.log(N=N, TT=TT, y=data, W=W)
ols_log_2chi

estim_log2_chi <- log_lin_estimnar2(x0=ols_log_2chi, N=N, TT=TT, y=data, W=W)
estim_log2_chi


round(estim_log1_chi$score,4)    # estimation problems 
sum(estim_log1_chi$coeflog[-1])  # estimates at the boundary of parameter space

## try to use same optimizer without stationarity constraint

optz <- list("algorithm"="NLOPT_LD_SLSQP", "xtol_rel"=1.0e-8)
no_staz_log1 <- nloptr(x0=ols_log_chi, eval_f=logl_log_lin, eval_grad_f=scor_log,
                  #eval_g_ineq=constr, eval_jac_g_ineq=j_constr,
                  opts=optz, N=N, TT=TT, y=data, W=W)

no_staz_log1
scor_log(no_staz_log1$solution, N, TT, data, W) # <- works fine 
sum(no_staz_log1$solution[-1])                  # <- but parameters out of stationary region

no_staz_log2 <- nloptr(x0=ols_log_2chi, eval_f=logl_log_lin2, eval_grad_f=scor_log2,
                  #eval_g_ineq=constr, eval_jac_g_ineq=j_constr,
                  opts=optz, N=N, TT=TT, y=data, W=W)

no_staz_log2
scor_log2(no_staz_log2$solution, N, TT, data, W) # <- works ok 
sum(no_staz_log2$solution[-1])                   # <- but again out of stationary region


## testing linearity

####### Testing linearity versus ID model
estim_lin1_chi

ID_test <- score_test_nonlin_h0(c=estim_lin1_chi$coeflin, N=N, TT=TT, y=data, W=W)
ID_test 
qchisq(p=0.10, df=1, lower.tail = FALSE)
qchisq(p=0.05, df=1, lower.tail = FALSE)
qchisq(p=0.01, df=1, lower.tail = FALSE)

# Result: linearity rejected at all usual significance levels



########### Testing linearity versus STNAR model

# optimize test with respect to nuisance parameter gamma
optg_stnar <- optim_gamma_stpnar(par=estim_lin1_chi$coeflin, N=N, TT=TT, y=data, W=W)
optg_stnar

# take extremes of grid of values in the neighborhood of the optimum
gam_L1 <- max(0.0001, optg_stnar$minimum-0.01)
gam_U1 <- optg_stnar$minimum+0.01

coef.lin.stnar <- c(estim_lin1_chi$coeflin, 0, 0) # <- vector of estimated parameters under the null hypothesis
LMstnar_opt <- score_test_star1_j(coef.lin.stnar, N, TT, y=data, W, J=99, l=10, gam_L=gam_L1, gam_U=gam_U1)
LMstnar_opt

# Result:


# select arbitrarily extremes of grid of values
gam_L <- 0.005
gam_U <- 3

LMstnar_a <- score_test_star1_j(coef.lin.stnar, N, TT, y=data, W, J=99, l=10, gam_L=gam_L, gam_U=gam_U)
LMstnar_a

# Result:

############ Testing linearity versus TNAR model

# optimize test with respect to nuisance parameter gamma
optg_tnar <- optim_gamma_tpnar(par=estim_lin1_chi$coeflin, N=N, TT=TT, y=data, W=W)
optg_tnar

# take extremes of grid of values in the neighborhood of the optimum
Tgam_L2 <- max(0.0001, optg_tnar$minimum-0.01)
Tgam_U2 <- optg_tnar$minimum+0.01

coef.lin.tnar <- c(estim_lin1_chi$coeflin, 0, 0, 0, 0) # <- vector of estimated parameters under the null hypothesis
LMtnar_opt <- score_test_tar1_j(coef.lin.tnar, N, TT, y=data, W, J=1, l=10, gam_L=Tgam_L2, gam_U=Tgam_U2)
LMtnar_opt

# Results:


# select arbitrarily extremes of grid of values

X_data=W%*%data
qqx <- matrix(0, 2, nrow(X_data))
for (k in 1:nrow(X_data)){
  qqx[,k] <- as.matrix(quantile(X_data[k,], prob=c(0.10, 0.90)))
}
Tgam_L <- max(0.0001, min(qqx[1,]))
Tgam_U <- max(qqx[2,])

LMtnar_a <- score_test_tar1_j(coef.lin.tnar, N, TT, y=data, W, J=1, l=10, gam_L=Tgam_L, gam_U=Tgam_U)
LMtnar_a

# Result: 