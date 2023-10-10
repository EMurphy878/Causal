setwd("~/Desktop/EDMS647 - Causal")
library(ggplot2)
dat<-read.csv("nclb.csv")
dim(dat)
names(dat)
str(dat)
head(dat) 
tail(dat)
#Proficiency Categories 
dat$prof3 <- cut(dat$s.prof, c(0, 50, 75, 100), labels = c('low', 'medium', 'high'), 
                 right = F)
#Standard Categories 
dat$standard <- 100 - dat$s.prof
# standards categories
dat$stand3 <- cut(dat$standard, c(0, 25, 50, 100), labels = c('low', 'medium', 'high'))
table(dat$prof3, dat$stand3)

#Date of Policy Intervention (From 2003 after)
dat$policy <- ifelse(dat$year > 2002, 1, 0)
#Centered Time Variable at 2002
dat$year02 <- dat$year - 2002

dat <- subset(dat, state != 'New York' & state != 'Vermont')

#################Part 2 ##################
### Model 1###
mod.1 <-as.formula(math8~ policy*year02 )
out.1 <- lm(mod.1, data=dat)
summary(out.1)

#Predict NCLB Effect for 2011 (year02 = 9)
predict(out.1, data.frame(policy = 1, year02 = 9)) - 
  predict(out.1, data.frame(policy = 0, year02 = 9))

######Plot - FINISH#####
m.tab <- with(dat, tapply(math8, year, mean, na.rm = TRUE))
m.tab <- m.tab[!is.na(m.tab)]   # remove NAs
NewDF <- data.frame(cbind(Year = as.numeric(as.vector(names(m.tab))),
                          Math8 = as.numeric(m.tab)))
plot(Math8~Year, data = NewDF, pch = 20,  xlab = "Year", ylab = "Math 8 Score",
     main = 'Math Achievement Scores (8th Graders)')
abline(v = 2002)
lines(as.numeric(names(m.tab[1:4])), m.tab[1:4], type = 'o', lwd = 2)
lines(as.numeric(names(m.tab[5:9])), m.tab[5:9], type = 'o', lwd = 2)







######### Part 3 ##########
mod.2 <- as.formula(math8~ stand3*year02*policy )
out.2 <-lm(mod.2, data = dat)
summary(out.2)


#Estimate Effect for 2011

################## Part 4#########
#State Fixed Effects#
mod.3 <-(math8 ~ year02:stand3 + policy:stand3 + policy:stand3:year02 + state)
out.3 <-lm(mod.3, data = dat)
summary(out.3)

#Plot
m.tab <- with(dat, tapply(math8, list(year,stand3), mean, na.rm = TRUE))

plot(math8 ~ year, data = dat, cex = .7, col = as.numeric(stand3),
     main = 'Math Achievement Scores (8th Graders)')
abline(v = 2002)
for (i in 1:3) {
  m.tab <- with(dat[dat$stand3 == levels(dat$stand3)[i], ], tapply(math8, year, mean, na.rm = TRUE))
  m.tab <- m.tab[!is.na(m.tab)]   # remove NAs
  lines(as.numeric(names(m.tab[1:4])), m.tab[1:4], col = i, type = 'o', lwd = 2)
  lines(as.numeric(names(m.tab[5:9])), m.tab[5:9], col = i, type = 'o', lwd = 2)
}
#2011 Prediction
predict(out.3, data.frame(policy = 1, year02 = 9)) - 
  predict(out.3, data.frame(policy = 0, year02 = 9))

#State and Time Fixed Effects#
mod.4 <-(math8 ~ year02:stand3 + policy:stand3 + policy:stand3:year02 + state + year02)
out.4 <-lm(mod.4, data = dat)
summary(out.4)
