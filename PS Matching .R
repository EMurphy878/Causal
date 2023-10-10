setwd("~/Desktop/EDMS647 - Causal")
library(sem)
library(ggplot2)
dat <- read.csv("analysis1_dataset.csv")

################################# Part 1 - Assess Initial Imbalance############################# 
#Descriptive statistics by treatment group 
#Mean Difference between Treatment and Control 
  #Note: May need to remove 'other','age','married','momdegr','daddegr','credit'
  #Note vm = treatment variable,  0 = vocab (control), 1 = math (treatment)
#All Baseline Covariates 
vars <- c('mathpre', 'vocabpre', 'collgpaa', 'hsgpaar', 'actcomp', 'credit', 
          'majormi', 'male', 'age', 'afram', 'cauc', 'beck', 'mars', 'pintell', 
          'pemot', 'pconsc', 'pagree', 'pextra', 'preflit', 'likelit', 'likemath',
          'numbmath','other','age','married','momdegr','daddegr','credit')

dat$vm.f <- ifelse(dat$vm == '0', 1, 0) 

table(dat$vm.f)
dat$z <- factor(dat$vm)
relevel(dat$z, ref = "1")
table(dat$z)

#Check Initial Imbalance 
#Add Standard Mean Difference (SMD using Glass's Delta method) and Variance Ratio (VR) to Balance 
select.diff <- function(x, grp, wt = rep(1, length(x))) 
{
  if (is.factor(grp)) grp <- as.numeric(grp) - 1     # transform dichotomeous factor into dummy
  if (is.factor(x)) x <- as.numeric(x) - 1           # transform dichotomeous factor into dummy
  stat.wt <- function(xw) {                          # function computing weighted mean and variance
    out <- cov.wt(xw[, 1, drop = F], xw[, 2])
    c(mean.wt = out$center, var.wt = out$cov)
  }
  stat <- unlist(by(cbind(x, wt), grp, stat.wt))     # weighted mean and variance (vector)
  m <- stat[seq(1, 4, by = 2)]                       # extract means
  v <- stat[seq(2, 4, by = 2)]                       # extract variances
  B <- (m[2] - m[1]) / sqrt(sum(v) / 2)              # standardized difference in means
  R <- v[2] / v[1]                                   # variance ratio
  rslt <- summary(lm(x ~ grp, weights = wt))         # run regression test
  res <- c(rslt$coef[2, c(1:2, 4)], B, R)
  names(res) <- c("Mean.Diff", "Std.Error", "p-value", "Std.Mean.Diff", "Var.Ratio")
  res
}

#Initial Imbalance in Covariates 
imbal1 <- t(sapply(dat[, vars], select.diff, grp = dat$vm.f))
imbal1
B <- imbal1[, "Std.Mean.Diff"]
R <- imbal1[, "Var.Ratio"]
plot(B, R, xlim = c(-1, 1), pch = 16,
     main = 'Imbalance', xlab = 'Std. Mean Difference', ylab = 'Variance Ratio')
abline(h = 1, v = 0)
abline(h = c(4/5, 5/4), v = c(-.1, .1), lty = 3)
ind <- (abs(B) > .1 | R < 4/5 | R > 5/4) 
text(B[ind], R[ind], vars[ind], pos = 2 + 2*(B[ind] > 0), offset = .4, cex = .7, xpd = T)

#Kernel - Density Plots for each variable separated by treatment group 
dnsty <- function(x, grp, wt = rep(1, length(x)), main = NULL, 
                  lwd = 1, col.vec = c('blue', 'red'), ...) 
{
  # plots kernel-density estimates for groups
  # x   ... metric variable
  # grp ... grouping variable (factor)
  # wt  ... weights for observations (e.g., inverse-propensity weights or survey weights)
  
  lev <- levels(grp)
  cmplt <- na.omit(data.frame(x, grp, wt))
  x <- cmplt$x
  grp <- cmplt$grp
  wt <- cmplt$wt
  dens <- list()
  x.m <- rep(NA, length(table(grp)))
  bw <- density(x, weights = wt / sum(wt))$bw
  for (i in 1:length(table(grp))) {
    grp.ind <- grp == lev[i]
    dens[[i]] <- if (is.null(wt)) density(x[grp == lev[i]], ...) else
      density(x[grp.ind], na.rm = T, weights = wt[grp.ind] / sum(wt[grp.ind]), bw = bw, ...)
    x.m[i] <- weighted.mean(x[grp.ind], w = wt[grp.ind]) 
  }
  plot(dens[[1]], type = 'n', main = main, xlab = '', ylab = '', yaxt = 'n',
       ylim = c(0, max(c(dens[[1]]$y, dens[[2]]$y))), 
       xlim = range(c(dens[[1]]$x, dens[[2]]$x)), bty = 'l')
  for (i in 1:length(table(grp))) {
    grp.ind <- grp == lev[i]
    ind <- dens[[i]]$x - x.m[i]
    lines(dens[[i]], col = col.vec[i], lwd = lwd, lty = 3 - i, xpd = T)
    lines(rep(x.m[i], 2), c(0, dens[[i]]$y[abs(ind) == min(abs(ind))][1]), col = col.vec[i], lwd = lwd, lty = 3)
  }
}

for (i in 1:length(vars))
{
  if ((i-1) %% 6 == 0) {par(mfrow = c(2, 3))} # creates a plotting window for 6 plots
  nam <- vars[i]
  x <- dat[, nam]
  if (!is.numeric(x)) x <- as.numeric(x)
  dnsty(x, dat$z, main = nam)
}


##########Check Overlap 
overlap <- function(x, z, lab = NULL, bin = 20)
{
  # plot a histogram of a covariate by group
  # x   ... numeric vector (covariate)
  # z   ... treatment indicator (dummy)
  # lab ... label for title and x-axis
  # bin ... number of bins for histogram
  
  r1 <- range(x)
  if (!is.numeric(z)) z <- as.numeric(z) - 1
  c.dat <- hist(x[z == 0], seq(r1[1], r1[2], length = bin), plot = F)  # histogram data for control group
  t.dat <- hist(x[z == 1], seq(r1[1], r1[2], length = bin), plot = F)  # histogram data for treatm. group
  t.dat$counts <- -t.dat$counts
  plot(c.dat, axes = F, ylim = c(min(t.dat$counts), max(c.dat$counts)),
       main = lab, xlab = lab)
  plot(t.dat, add = T, density = 30)
  axis(1)
  ax.loc <- axis(2, labels = F)
  axis(2, at = ax.loc, labels = abs(ax.loc))
  y <- par('usr')[3:4]
  text(rep(max(x), 2), c(y[2] - diff(y)*.05, y[1] + diff(y)*.05), 
       c('Control', 'Treatment'), adj = 1, xpd = T)
}

for (i in 1:length(vars))
{
  if((i-1) %% 8 == 0) {par(mfrow = c(2, 4))} # creates a plotting window for 6 plots
  nam <- vars[i]
  x <- dat[, nam]
  if (!is.numeric(x)) x <- as.numeric(x)
  overlap(x, dat$vm.f, nam)
}

######################## Part 2 - PS Estimation of ATT#############
#A - Inverse Propensity Score Weighting#

#Scatterplot Matrix
library(car)
scatterplotMatrix(~ logit.mod1 + mathpre + vocabpre, data = dat)

#Main Effects Model 
adj.set <- c("mathpre", "vocabpre", "actcomp", "hsgpaar", "collgpaa",
             "numbmath", "likemath", "likelit", "preflit", "majormi", "mars", 
             "cauc", "afram", "other", "male")
mod.1 <- as.formula(vm.f ~ mathpre + vocabpre + actcomp + hsgpaar + collgpaa +
                     numbmath + likemath + likelit + preflit + majormi +
                     mars + cauc + afram + male) 
out1 <- glm(mod.1, dat, family = binomial(link = logit))
dat$logit.mod1 <- log((out1$fitted.values)/(1-out1$fitted.values))

scatterplotMatrix(~logit.mod1 + mathpre +vocabpre, data = dat)
summary(out1)

#Model w/ Interaction Terms 
mod.3 <- as.formula(vm.f ~ mathpre + vocabpre + actcomp + hsgpaar*collgpaa +
                      numbmath*likemath + likelit*preflit+ majormi*numbmath +
                      mars + cauc + afram + male) 
out3 <- glm(mod.3, dat, family = binomial(link = logit))

summary(out3)

#PS With Mod 1 (Main Effects)
dat$ps1 <- out1$fitted.values
dat$logps1 <- log(dat$ps1 / (1-dat$ps1))
hist(dat$ps1, breaks = 30)
hist(dat$logps1, breaks = 30)
overlap(dat$ps1, dat$vm.f, bin = 30)
overlap(dat$logps1, dat$vm.f, bin = 30)

#PS with Mod 2 
dat$ps2 <- out2$fitted.values
dat$logps2 <- log(dat$ps2 / (1-dat$ps2))
hist(dat$ps2, breaks = 30)
hist(dat$logps2, breaks = 30)
overlap(dat$ps2, dat$vm.f, bin = 30)
overlap(dat$logps2, dat$vm.f, bin = 30)

#PS With Mod 3 
dat$ps3 <- out3$fitted.values
dat$logps3 <- log(dat$ps3 / (1-dat$ps3))

hist(dat$ps3, breaks = 30)
hist(dat$logps3, breaks = 30)
overlap(dat$ps3, dat$vm.f, bin = 30)
overlap(dat$logps3, dat$vm.f, bin = 30)

discard <- function(lps, grp, caliper = .05)
{
  # creates an index for discarding non-overlapping cases (with caliper)
  # returns logical vector (F = overlapping, T = non-overlapping)
  # lps     ... propensity score logit
  # grp     ... treatment indicator
  # caliper ... caliper in SD of lps
  
  ovl <- max(tapply(lps, grp, min)) - sd(lps) * caliper   # left limit of overlap
  ovr <- min(tapply(lps, grp, max)) + sd(lps) * caliper   # right limit of overlap
  ind <- lps > ovr | lps < ovl
  cat('Number of non-overlapping cases:', sum(ind), '::: breaks:', ovl, ovr, '\n')
  ind
}

(del.ind <- discard(dat$logps3, dat$vm.f, caliper = .05))

#Weights for computing ATT T= 1, C = PS/(1-PS) with Mod 1  
dat$ipwc1 <- with(dat, ps1/(1-ps1))
dat$ipw1 <- 1 
dat$ipw1 <- ifelse(dat$vm.f == 1,1, dat$ipwc1)
#Set weights of non overlapping cases to 0 
dat$ipwo1 <-dat$ipw1
dat$ipwo1[del.ind] <- 0


#Weights for computing ATT T= 1, C = PS/(1-PS) with Mod 3  
dat$ipwc3 <- with(dat, ps3/(1-ps3))
dat$ipw3 <- 1 
dat$ipw3 <- ifelse(dat$vm.f == 1,1, dat$ipwc3)
dat$ipwo3 <-dat$ipw3
dat$ipwo3[del.ind] <- 0

#Check Imbalance using Model 1 Weights Including Non-overlapping Cases *Use This One
imbal.mod1 <- t(sapply(dat[, vars], select.diff, grp = dat$vm.f, wt = dat$ipw1))
imbal.mod1
#Check Imbalance using Model 3 Weights Including Non-overlapping Cases 
imbal.mod3 <- t(sapply(dat[, vars], select.diff, grp = dat$vm.f, wt = dat$ipw3))

imbal.mod3
#Check Imbalance using Model 3 After Removing Non-overlapping Cases *Or this one
imbal.mod3.2 <- t(sapply(dat[, vars], select.diff, grp = dat$vm.f, wt = dat$ipwo3))
imbal.mod3.2
imbal1

#Plot Imbalance Check *Use Mod 3 
B <- imbal.mod3[, "Std.Mean.Diff"]
R <- imbal.mod3[, "Var.Ratio"]
plot(B, R, xlim = c(-1, 1), pch = 16,
     main = 'Imbalance Mod3', xlab = 'Std. Mean Difference', ylab = 'Variance Ratio')
abline(h = 1, v = 0)
abline(h = c(4/5, 5/4), v = c(-.1, .1), lty = 3)
ind <- (abs(B) > .1 | R < 4/5 | R > 5/4) 
text(B[ind], R[ind], vars[ind], pos = 2 + 2*(B[ind] > 0), offset = .4, cex = .7, xpd = T)

dnsty(dat$ps1, dat$vm.f, wt = dat$ipw1, main = "PS Logit Mod 1")
dnsty(dat$ps3, dat$z, wt = dat$ipw3, main = "PS Logit Mod 3")



for (i in 1:length(vars))
{
  if ((i-1) %% 6 == 0) {par(mfrow = c(2, 3))} # creates a plotting window for 6 plots
  nam <- vars[i]
  x <- dat[, nam]
  if (!is.numeric(x)) x <- as.numeric(x)
  dnsty(x, dat$vm.f, wt = dat$ipw3, main = nam)
}



#Estimate ATT Using Inverse Propensity Weighting 

#Biased Treatment Effect (without PS Adjustement)
summary(lm(mathall ~ vm.f, data = dat))
#Treatment Effect with PS ATT Adjustment (All Cases) Using Model 3 Weights
summary(lm(mathall ~ vm.f, data = dat, weights = ipw3))

#Treatment Effect with PS ATT Adjustment (w only Overlapping Cases) Using Model 3 Weights 
summary(lm(mathall ~ vm.f, data = dat, weights = ipwo3))

#Doubly Robust Estimation  *Using Model 3 Covariates and Weights 
mdl.3 <- as.formula(mathall ~ vm.f + mathpre + vocabpre + actcomp + hsgpaar + collgpaa +
                      numbmath + likemath + likelit + preflit+ majormi + numbmath +
                      mars + cauc + afram + male) 
summary(lm(mdl.3, data = dat, weights = ipw3))
summary(lm(mdl.3, data = dat, weights = ipwo3)) 


########### Part 3 - Propensity Score Stratification ########################

#Stratum for PS Stratification by Quantiles using PS3 Weights  
dat$ps5 <- cut(dat$ps3, quantile(dat$ps3, probs = seq(0,1, by = .2)), include.lowest = TRUE)
which(is.na(dat$ps5));
overlap(dat$ps3, dat$vm.f, 30)
abline(v = quantile(dat$ps3, seq(0, 1, by = .2)), col = 'red', lty = 2)
#Compute ATT Weights 
#USE THIS FOR Stratum Weights 
dat$z <- factor(dat$vm)
relevel(dat$z, ref = "1")
table(dat$z)
ps.strat <- table(dat$z, dat$ps5)
O <- table(dat$z, dat$ps5)
E <- rbind(ps.strat[1,],131*ps.strat[1,]/79)
W <- E/O
dat$strwt<- W[cbind(dat$z,dat$ps5)] 


O <- table(dat$z, dat$ps5)
E <- rbind(O[1,],131*O[1,]/79)
W <- E/O
dat$strwt<- W[cbind(dat$z,dat$ps5)] 


#Check Balance 
#Variance Ratio and Std Mean Diff 
imbal <- t(sapply(dat[, vars], select.diff, grp = dat$vm.f, wt = dat$strwt))

B <- imbal[, "Std.Mean.Diff"]
R <- imbal[, "Var.Ratio"]
plot(B, R, xlim = c(-1, 1), pch = 16,
     main = 'Imbalance', xlab = 'Std. Mean Difference', ylab = 'Variance Ratio')
abline(h = 1, v = 0)
abline(h = c(4/5, 5/4), v = c(-.1, .1), lty = 3)
ind <- (abs(B) > .1 | R < 4/5 | R > 5/4) 
text(B[ind], R[ind], vars[ind], pos = 2 + 2*(B[ind] > 0), offset = .4, cex = .7, xpd = T)

#Kernel Density Plot
for (i in 1:length(vars))
{
  if ((i-1) %% 6 == 0) {par(mfrow = c(2, 3))} # creates a plotting window for 6 plots
  nam <- vars[i]
  x <- dat[, nam]
  if (!is.numeric(x)) x <- as.numeric(x)
  dnsty(x, dat$vm.f, wt = dat$strwt, main = nam)
}

#Estimate Treatment Effect  
summary(lm(mathall ~ vm.f, data = dat, weights = strwt))
#Doubly Robust Treatment Effect 
summary(lm(mdl.3, data =dat, weights = strwt))

#################### Part 4 - Matching####################
install.packages('optmatch')
library(optmatch)

#ID Matches 
ps.dist <- match_on(out3)
dim(ps.dist)                                   
table(dat$vm.f) 
ps.dist[1:10, 1:5]  
(dat$mtch <- fullmatch(ps.dist, data = dat)) 
summary(dat$mtch)                              # summary of matching process
table(dat$mtch)                                # each indicator refers to one treatment and one control case
length(table(dat$mtch))                        # number of matched strata 
table(dat$mtch, scs$z)

#Weight Matches 
(O <- table(dat$vm.f, dat$mtch))                                        # observed table
(E <- outer(apply(O, 1, sum), O['1', ] / sum(O['1', ])))                # expected table
W <- E / O                                                        # compute weights
dat$fmwt <- W[cbind(as.character(dat$vm.f), as.character(dat$mtch))]  

dat$ipwo3[del.ind] <- 0
dat$fmwt2 <- dat$fmwt
dat$fmwt2 [del.ind] <- 0
#Check Balance 
imbal <- t(sapply(dat[, vars], select.diff, grp = dat$vm.f, wt = dat$fmwt))

B <- imbal[, "Std.Mean.Diff"]
R <- imbal[, "Var.Ratio"]
plot(B, R, xlim = c(-1, 1), pch = 16,
     main = 'Imbalance', xlab = 'Std. Mean Difference', ylab = 'Variance Ratio')
abline(h = 1, v = 0)
abline(h = c(4/5, 5/4), v = c(-.1, .1), lty = 3)
ind <- (abs(B) > .1 | R < 4/5 | R > 5/4) 
text(B[ind], R[ind], vars[ind], pos = 2 + 2*(B[ind] > 0), offset = .4, cex = .7, xpd = T)


#Kernel Density Plot
for (i in 1:length(vars))
{
  if ((i-1) %% 6 == 0) {par(mfrow = c(2, 3))} # creates a plotting window for 6 plots
  nam <- vars[i]
  x <- dat[, nam]
  if (!is.numeric(x)) x <- as.numeric(x)
  dnsty(x, dat$z, wt = dat$fmwt, main = nam)
}

#Estimate Treatment Effect 
summary(lm(mathall ~ vm.f, data = dat, weights = fmwt))
#Doubly Robust Treatment Effect 
summary(lm(mdl.3, data =dat, weights = fmwt))

###### Part 5 - Regression Estimates ########################
balance.reg <- function(z, y, lps)
{
  # extract only required stats from regression test
  # z   ... independent variable (treatment indicator)
  # y   ... dependent variable (covariates)
  # lps ... PS-logit
  
  if (is.factor(y)) y <- as.numeric(y) - 1       # transform dichotomeous factor into dummy
  rslt <- summary(lm(y ~ z + lps + I(lps^2) + I(lps^3)))$coef[2, c(1:2, 4)]  # run regression test
  grp.var <- tapply(y, z, var)                   # group variances (unadjusted!)
  rslt <- c(rslt, rslt[1] / sqrt(mean(grp.var)), grp.var[2] / grp.var[1])
  names(rslt) <- c("Mean Diff.", "Std. Error", "p-value", "Cohen's d", "Var. Ratio")
  round(rslt, 4)
}

# apply balance.reg() to all covariates: "balance statistics"
# note: variance ratios are not adjusted for the regression adjustment!
(imbal <- sapply(dat[, vars], balance.reg, z = dat$vm.f, lps = dat$logps3))

# ::::: plot standardized mean difference & variance ratio :::::
plot(imbal["Cohen's d", ], imbal["Var. Ratio", ], xlim = c(-1, 1),
     main = 'Imbalance', xlab = 'Std. Mean Difference', ylab = 'Variance Ratio')
abline(h = 1, v = 0)
abline(h = c(4/5, 5/4), v = c(-.1, .1), lty = 3)

#Estimate of ATT 
out.v <- lm(mathall ~ logps3 + I(logps3^2) + I(logps3^3), data = dat, subset = vm.f == 0)
summary(out.m)
pot.v <- predict(out.v, dat)
dat$cdiff <- dat$mathall-pot.v
sum(dat$cdiff)/79

#Std Error for ATT Regression Estimate
library(boot)
PSreg <- function(bdat, i)
{
  # draw bootstrap sample (done by the boot() function)
  bsmpl <- bdat[i, ] 
  # PS model
  out1 <- glm(mod.3, data = bsmpl, family = 'binomial')
  # get PS, PS-logit 
  bsmpl$ps <- out1$fitted                      # fitted values are the PS
  bsmpl$lps <- log(bsmpl$ps / (1 - bsmpl$ps))  # PS-logit = log(PS/(1-PS))
  out.v <- lm(mathall ~ lps + I(lps^2) + I(lps^3), data = bsmpl, subset = vm.f == 0)
  pot.v <- predict(out.v, bsmpl)   # predicted "potential" control outcomes
  rslt <- sum(dat$mathall-pot.v)/79      # treatment effect
  rslt
}

(bootstat <- boot(dat, statistic = PSreg, R = 1000, stype = 'i'))

#Doubly Roubust Estimte 
mdl3 <- as.formula(mathall ~ mathpre + vocabpre + actcomp + hsgpaar + collgpaa +
                     numbmath + likemath + likelit + preflit + majormi +
                     mars + cauc + afram + male +
                     I(logps3^2) + I(logps3^3))

out.v2 <- lm(mdl3, data = dat, subset = vm.f == 0)
pot.v2 <- predict(out.v2, dat)
dat$cdiff2 <- dat$mathall-pot.v2
sum(dat$cdiff2)/79

PSreg <- function(bdat, i)
{
  # draw bootstrap sample (done by the boot() function)
  bsmpl <- bdat[i, ] 
  # PS model
  out1 <- glm(mod.3, data = bsmpl, family = 'binomial')
  # get PS, PS-logit 
  bsmpl$ps <- out1$fitted                      # fitted values are the PS
  bsmpl$lps <- log(bsmpl$ps / (1 - bsmpl$ps))  # PS-logit = log(PS/(1-PS))
  out.v <- lm(mdl3, data = bsmpl, subset = vm.f == 0)
  pot.v <- predict(out.v, bsmpl)   # predicted "potential" control outcomes
  rslt <- sum(dat$mathall-pot.v)/79      # treatment effect
  rslt
}

(bootstat <- boot(dat, statistic = PSreg, R = 1000, stype = 'i'))

library(CBPS)
out.cb <- CBPS(mod.3, data = dat, ATT = 1)
balance(out.cb)
plot(out.cb)
