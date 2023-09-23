
###########################################################################
## Functions to visualize simulations results of
## "SimulationStudyPNAR.Rdata" or "SimulationStudyLogPNAR.Rdata"
###########################################################################

### Visualize mean of estimated parameters over the simulations
### with related standard errors

#######Small N, large TT
round(colMeans(out1),3)
round(colMeans(out2),3)
round(colMeans(out3),3)
round(colMeans(out4),3)
round(colMeans(out5),3)
round(colMeans(out6),3)
#######Large N, large TT
round(colMeans(out7),3)
round(colMeans(out8),3)
round(colMeans(out9),3)
round(colMeans(out10),3)
round(colMeans(out11),3)
round(colMeans(out12),3)
#######Large N, small TT
round(colMeans(out13),3)
round(colMeans(out14),3)
round(colMeans(out15),3)
round(colMeans(out16),3)
round(colMeans(out17),3)
round(colMeans(out17),3)

############################################################
############################################################
##  freq() --- function computing the frequency of cases
##  where the coefficient is significant at 5% significance
##  level, according to z-test
##  Inputs:
##    x = values to test
##    S = number of replications
##  output:
##    frequency of z-tests >1.96 over S replications
############################################################

freq<- function(x,S){
  sum(x>1.96)/S*100
}

####### frequence of cases the coefficient is significant at 5%

apply(out1[,7:9], 2, FUN=freq, S=1000)
apply(out1[,26:30], 2, FUN=freq, S=1000)

apply(out2[,7:9], 2, FUN=freq, S=1000)
apply(out2[,26:30], 2, FUN=freq, S=1000)

apply(out3[,7:9], 2, FUN=freq, S=1000)
apply(out3[,26:30], 2, FUN=freq, S=1000)

apply(out4[,7:9], 2, FUN=freq, S=1000)
apply(out4[,26:30], 2, FUN=freq, S=1000)

apply(out5[,7:9], 2, FUN=freq, S=1000)
apply(out5[,26:30], 2, FUN=freq, S=1000)

apply(out6[,7:9], 2, FUN=freq, S=1000)
apply(out6[,26:30], 2, FUN=freq, S=1000)

apply(out7[,7:9], 2, FUN=freq, S=1000)
apply(out7[,26:30], 2, FUN=freq, S=1000)

apply(out8[,7:9], 2, FUN=freq, S=1000)
apply(out8[,26:30], 2, FUN=freq, S=1000)

apply(out9[,7:9], 2, FUN=freq, S=1000)
apply(out9[,26:30], 2, FUN=freq, S=1000)

apply(out10[,7:9], 2, FUN=freq, S=1000)
apply(out10[,26:30], 2, FUN=freq, S=1000)

apply(out11[,7:9], 2, FUN=freq, S=1000)
apply(out11[,26:30], 2, FUN=freq, S=1000)

apply(out12[,7:9], 2, FUN=freq, S=1000)
apply(out12[,26:30], 2, FUN=freq, S=1000)

apply(out13[,7:9], 2, FUN=freq, S=1000)
apply(out13[,26:30], 2, FUN=freq, S=1000)

apply(out14[,7:9], 2, FUN=freq, S=1000)
apply(out14[,26:30], 2, FUN=freq, S=1000)

apply(out15[,7:9], 2, FUN=freq, S=1000)
apply(out15[,26:30], 2, FUN=freq, S=1000)

apply(out16[,7:9], 2, FUN=freq, S=1000)
apply(out16[,26:30], 2, FUN=freq, S=1000)

apply(out17[,7:9], 2, FUN=freq, S=1000)
apply(out17[,26:30], 2, FUN=freq, S=1000)

apply(out18[,7:9], 2, FUN=freq, S=1000)
apply(out18[,26:30], 2, FUN=freq, S=1000)


############# frequency of cases where information criteria select the right model

x1.a <- out1[,c(13,36)]==apply(out1[,c(13,36)], 1, min)
100*apply(x1.a,2,mean)

x1.b <- out1[,c(14,37)]==apply(out1[,c(14,37)], 1, min)
100*apply(x1.b,2,mean)

x1.q <- out1[,c(15,38)]==apply(out1[,c(15,38)], 1, min)
100*apply(x1.q,2,mean)


x2.a <- out2[,c(13,36)]==apply(out2[,c(13,36)], 1, min)
100*apply(x2.a,2,mean)

x2.b <- out2[,c(14,37)]==apply(out2[,c(14,37)], 1, min)
100*apply(x2.b,2,mean)

x2.q <- out2[,c(15,38)]==apply(out2[,c(15,38)], 1, min)
100*apply(x2.q,2,mean)


x3.a <- out3[,c(13,36)]==apply(out3[,c(13,36)], 1, min)
100*apply(x3.a,2,mean)

x3.b <- out3[,c(14,37)]==apply(out3[,c(14,37)], 1, min)
100*apply(x3.b,2,mean)

x3.q <- out3[,c(15,38)]==apply(out3[,c(15,38)], 1, min)
100*apply(x3.q,2,mean)


x4.a <- out4[,c(13,36)]==apply(out4[,c(13,36)], 1, min)
100*apply(x4.a,2,mean)

x4.b <- out4[,c(14,37)]==apply(out4[,c(14,37)], 1, min)
100*apply(x4.b,2,mean)

x4.q <- out4[,c(15,38)]==apply(out4[,c(15,38)], 1, min)
100*apply(x4.q,2,mean)


x5.a <- out5[,c(13,36)]==apply(out5[,c(13,36)], 1, min)
100*apply(x5.a,2,mean)

x5.b <- out5[,c(14,37)]==apply(out5[,c(14,37)], 1, min)
100*apply(x5.b,2,mean)

x5.q <- out5[,c(15,38)]==apply(out5[,c(15,38)], 1, min)
100*apply(x5.q,2,mean)


x6.a <- out6[,c(13,36)]==apply(out6[,c(13,36)], 1, min)
100*apply(x6.a,2,mean)

x6.b <- out6[,c(14,37)]==apply(out6[,c(14,37)], 1, min)
100*apply(x6.b,2,mean)

x6.q <- out6[,c(15,38)]==apply(out6[,c(15,38)], 1, min)
100*apply(x6.q,2,mean)


x7.a <- out7[,c(13,36)]==apply(out7[,c(13,36)], 1, min)
100*apply(x7.a,2,mean)

x7.b <- out7[,c(14,37)]==apply(out7[,c(14,37)], 1, min)
100*apply(x7.b,2,mean)

x7.q <- out7[,c(15,38)]==apply(out7[,c(15,38)], 1, min)
100*apply(x7.q,2,mean)


x8.a <- out8[,c(13,36)]==apply(out8[,c(13,36)], 1, min)
100*apply(x8.a,2,mean)

x8.b <- out8[,c(14,37)]==apply(out8[,c(14,37)], 1, min)
100*apply(x8.b,2,mean)

x8.q <- out8[,c(15,38)]==apply(out8[,c(15,38)], 1, min)
100*apply(x8.q,2,mean)


x9.a <- out9[,c(13,36)]==apply(out9[,c(13,36)], 1, min)
100*apply(x9.a,2,mean)

x9.b <- out9[,c(14,37)]==apply(out9[,c(14,37)], 1, min)
100*apply(x9.b,2,mean)

x9.q <- out9[,c(15,38)]==apply(out9[,c(15,38)], 1, min)
100*apply(x9.q,2,mean)

x10.a <- out10[,c(13,36)]==apply(out10[,c(13,36)], 1, min)
100*apply(x10.a,2,mean)

x10.b <- out10[,c(14,37)]==apply(out10[,c(14,37)], 1, min)
100*apply(x10.b,2,mean)

x10.q <- out10[,c(15,38)]==apply(out10[,c(15,38)], 1, min)
100*apply(x10.q,2,mean)

x11.a <- out11[,c(13,36)]==apply(out11[,c(13,36)], 1, min)
100*apply(x11.a,2,mean)

x11.b <- out11[,c(14,37)]==apply(out11[,c(14,37)], 1, min)
100*apply(x11.b,2,mean)

x11.q <- out11[,c(15,38)]==apply(out11[,c(15,38)], 1, min)
100*apply(x11.q,2,mean)


x12.a <- out12[,c(13,36)]==apply(out12[,c(13,36)], 1, min)
100*apply(x12.a,2,mean)

x12.b <- out12[,c(14,37)]==apply(out12[,c(14,37)], 1, min)
100*apply(x12.b,2,mean)

x12.q <- out12[,c(15,38)]==apply(out12[,c(15,38)], 1, min)
100*apply(x12.q,2,mean)


x13.a <- out13[,c(13,36)]==apply(out13[,c(13,36)], 1, min)
100*apply(x13.a,2,mean)

x13.b <- out13[,c(14,37)]==apply(out13[,c(14,37)], 1, min)
100*apply(x13.b,2,mean)

x13.q <- out13[,c(15,38)]==apply(out13[,c(15,38)], 1, min)
100*apply(x13.q,2,mean)


x14.a <- out14[,c(13,36)]==apply(out14[,c(13,36)], 1, min)
100*apply(x14.a,2,mean)

x14.b <- out14[,c(14,37)]==apply(out14[,c(14,37)], 1, min)
100*apply(x14.b,2,mean)

x14.q <- out14[,c(15,38)]==apply(out14[,c(15,38)], 1, min)
100*apply(x14.q,2,mean)


x15.a <- out15[,c(13,36)]==apply(out15[,c(13,36)], 1, min)
100*apply(x15.a,2,mean)

x15.b <- out15[,c(14,37)]==apply(out15[,c(14,37)], 1, min)
100*apply(x15.b,2,mean)

x15.q <- out15[,c(15,38)]==apply(out15[,c(15,38)], 1, min)
100*apply(x15.q,2,mean)


x16.a <- out16[,c(13,36)]==apply(out16[,c(13,36)], 1, min)
100*apply(x16.a,2,mean)

x16.b <- out16[,c(14,37)]==apply(out16[,c(14,37)], 1, min)
100*apply(x16.b,2,mean)

x16.q <- out16[,c(15,38)]==apply(out16[,c(15,38)], 1, min)
100*apply(x16.q,2,mean)


x17.a <- out17[,c(13,36)]==apply(out17[,c(13,36)], 1, min)
100*apply(x17.a,2,mean)

x17.b <- out17[,c(14,37)]==apply(out17[,c(14,37)], 1, min)
100*apply(x17.b,2,mean)

x17.q <- out17[,c(15,38)]==apply(out17[,c(15,38)], 1, min)
100*apply(x17.q,2,mean)


x18.a <- out18[,c(13,36)]==apply(out18[,c(13,36)], 1, min)
100*apply(x18.a,2,mean)

x18.b <- out18[,c(14,37)]==apply(out18[,c(14,37)], 1, min)
100*apply(x18.b,2,mean)

x18.q <- out18[,c(15,38)]==apply(out18[,c(15,38)], 1, min)
100*apply(x18.q,2,mean)


############## q-q plots

# o107 <- (out1[,1]-mean(out1[,1]))/sd(out1[,1])
# o108 <- (out1[,2]-mean(out1[,2]))/sd(out1[,2])
# o109 <- (out1[,3]-mean(out1[,3]))/sd(out1[,3])
# o57 <- (out2[,1]-mean(out2[,1]))/sd(out2[,1])
# o58 <- (out2[,2]-mean(out2[,2]))/sd(out2[,2])
# o59 <- (out2[,3]-mean(out2[,3]))/sd(out2[,3])

# o107 <- (out13[,1]-mean(out13[,1]))/sd(out13[,1])
# o108 <- (out13[,2]-mean(out13[,2]))/sd(out13[,2])
# o109 <- (out13[,3]-mean(out13[,3]))/sd(out13[,3])
# o57 <- (out14[,1]-mean(out14[,1]))/sd(out14[,1])
# o58 <- (out14[,2]-mean(out14[,2]))/sd(out14[,2])
# o59 <- (out14[,3]-mean(out14[,3]))/sd(out14[,3])

o107 <- (out13[,1]-mean(out13[,1]))/sd(out13[,1])
o108 <- (out13[,2]-mean(out13[,2]))/sd(out13[,2])
o109 <- (out13[,3]-mean(out13[,3]))/sd(out13[,3])
o57 <- (out7[,1]-mean(out7[,1]))/sd(out7[,1])
o58 <- (out7[,2]-mean(out7[,2]))/sd(out7[,2])
o59 <- (out7[,3]-mean(out7[,3]))/sd(out7[,3])

#pdf(file="qq_lin_new.pdf", width=7, height=7)
pdf(file="qq_log_new.pdf", width=7, height=7)
par(mfrow=c(3,2))
title <- expression(atop(Normal~QQ-plot~"for"~beta[0]~","~"T"~"="~20))
qqnorm(o107, main=title)
qqline(o107)
title <- expression(atop(Normal~QQ-plot~"for"~beta[0]~","~"T"~"="~100))
qqnorm(o57, main=title)
qqline(o57)
title <- expression(atop(Normal~QQ-plot~"for"~beta[1]~","~"T"~"="~20))
qqnorm(o108, main=title)
qqline(o108)
title <- expression(atop(Normal~QQ-plot~"for"~beta[1]~","~"T"~"="~100))
qqnorm(o58, main=title)
qqline(o58)
title <- expression(atop(Normal~QQ-plot~"for"~beta[2]~","~"T"~"="~20))
qqnorm(o109, main=title)
qqline(o109)
title <- expression(atop(Normal~QQ-plot~"for"~beta[2]~","~"T"~"="~100))
qqnorm(o59, main=title)
qqline(o59)
dev.off()


