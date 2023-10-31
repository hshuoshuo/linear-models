#Shuo Han U09953590
#MA575 Linear Models Spring2022 
#05/03/2022

#Part I.A – Univariate Data Analysis – Mean Testing
#Null Hypothesis: 𝐻0: 𝜇 = 𝜇0 = 30 
#The averCrashes fatalities per 100,000 population is 30
S <- read.csv("accident-1.csv")

alpha <- 0.05
n <- dim(S)[1]

Y <- S[1:n, 2]
X1Old <- S[1:n, 3]
X2Old <- S[1:n, 4]
X3Old <- S[1:n, 5]

sx1Old  <- sd(X1Old)
sx2Old  <- sd(X2Old)
sx3Old  <- sd(X3Old)

x1barOld <- mean(X1Old)
x2barOld <- mean(X2Old)
x3barOld <- mean(X3Old)

SEOld <- sx1Old/sqrt(n)

tcrit <- qt(alpha/2, df = n-1, lower.tail=F)

#margin of error
eps <- tcrit * SEOld

#claimed value of the mean
mu0 <- 300

# Confidence Interval
Low <- x1barOld - eps
Low
Upper <- x1barOld + eps
Upper

# Test Statistics
tstat <- (x1barOld - mu0)/SEOld
tstat

#p-value
pval <- 2*pt(tstat, df=n-1, lower.tail = F)
pval

#Summary
metric_name <-c("CI.lower", "CI.upper", "claimed.mean","T.stat", "p-value","alpha")
metric_val <- c(Low, Upper, mu0, tstat, pval, alpha)
options(digits =7) 
Summary <- data.frame(metric_name, metric_val)
Summary


#Part I.B – Univariate Data Analysis – Standard Deviation Testing
#Null hypothesis: 𝐻0: 𝜎 = 𝜎0
alpha <- 0.05

#claimed value of the standard deviation
sd0 = 30

# Confidence Interval
LowC<-qchisq(alpha/2, df=n-1, lower.tail = T)
UpperC<-qchisq(alpha/2, df=n-1, lower.tail = F)
LowC
UpperC

# Test Statistics
tstatC<- (n-1)*(sx1Old/sd0)^2
tstatC

#p-value
pvalC <- 2*pt(abs(tstatC), df=n-1, lower.tail = F)
pvalC

#Summary
metric_nameC <-c("CI.lower","CI.upper", "claimed.sd","T.stat", "p-value","alpha")
metric_valC <- c(LowC, UpperC, sd0, tstatC, pvalC, alpha)

options(digits =7) 
SummaryC <- data.frame(metric_nameC, metric_valC)
SummaryC



#Part I.C – Normality Testing
Q1 <- qqnorm(X1Old, ylab = "Quantiles of Deaths", main = "NQQ plot of Deaths")
qqline(X1Old, col="orange", lwd=3)

cor(Q1$x,Q1$y)

S2 <- read.csv("accident-2.csv")

n2 <- dim(S2)[1]
X21Old <- S2[1:n2, 3]
sx21Old  <- sd(X21Old)
x21barOld <- mean(X21Old)

Q2 <- qqnorm(X21Old, ylab = "Quantiles of Deaths 2", main = "NQQ plot of Death 2")
qqline(X21Old, col="orange", lwd=3)

cor(Q2$x,Q2$y)

#Part I.D – Parameter Comparisons for Means
#Null hypothesis: 𝐻0: 𝜇1 = 𝜇2
#test statistic
xbarD<-x1barOld-x21barOld
SED<-sqrt((sx1Old^2/n)+(sx21Old/n2))

tcritD<-qt(alpha/2, n+n2-2,lower.tail=F)
epsD<-tcritD*SED
tstatD<-(xbarD-0)/SED

#p-value calculated
pvalD<-2*pt(-abs(tstatD), n+n2-2, lower.tail = T)

#Confidence interval
LowD<-xbarD-epsD
UpperD<-xbarD+epsD

#Summary
metric_nameD <-c("CI.lower", "CI.upper","T.stat","p-value", "alpha")
metric_valD<- c(LowD, UpperD, tstatD, pvalD, alpha)
options(digits =4) 
DataSummaryD <- data.frame(metric_nameD, metric_valD)
DataSummaryD


#Part I.E – Parameter Comparisons for Variances
source("nemolm2.r")
#Null hypothesis: 𝐻 : 𝜎2 = 𝜎2
#test statistic
sx21Old  <- sd(X21Old)
fstatV<-sx1Old^2/sx21Old^2

#Confidence interval
fcritLV<-qf(alpha/2, n-1, n2-1,lower.tail = T )
fcritUV<-qf(alpha/2, n-1, n2-1,lower.tail = F )

#p-value calculated
fstatLV<-min(fstatV, 1/fstatV)
fstatUV<-max(fstatV, 1/fstatV)
pvalFV<-pf(fstatLV, n-1, n2-1, lower.tail = T) + pf(fstatUV, n-1, n2-1, lower.tail = F)

#Summary
metric_nameFV <-c("CI.lower","CI.upper", "T.stat", "p-value","alpha")
metric_valFV <- c(fcritLV, fcritUV, fstatV, pvalFV, alpha)
options(digits =7) 
SummaryFV <- data.frame(metric_nameFV, metric_valFV)
SummaryFV



#Part II.A – Simple Linear Regression
#Standardize
X1 <- (X1Old-x1barOld)/sx1Old
X2 <- (X2Old-x2barOld)/sx2Old
X3 <- (X3Old-x3barOld)/sx3Old

sx1  <- 1
sx2  <- 1
sx3  <- 1

x1bar <- 0
x2bar <- 0
x3bar <- 0

ybar <- mean(Y)
sy <- sd(Y)
covs <- cov(X1, Y)
rs <- cor(X1, Y)

SE <- sx1/sqrt(n)

plot(X1, Y, xlab = "Crashes", ylab = "Deaths", main = "Deaths v. Crashes")

beta1hat <- rs*sy/sx1
beta0hat <- ybar - beta1hat*x1bar

SSE <-  sy^2*(n-1)*(1-rs^2)
SE.beta1hat <- (1/sx1)*sqrt(SSE/(n-1)*(n-2))
SE.beta0hat <- sqrt(SSE/(n-2))*sqrt(1/n + (x1bar)^2/(sx1^2*(n-1)))

metric.name <-c("covariance", "r value", "r^2 value", "beta1hat", "SE.beta1hat", "beta0hat",
                "SE.beta0hat", "SSE")
metric.val <- c(covs, rs, rs^2, beta1hat, SE.beta1hat, beta0hat, SE.beta0hat, SSE)

D <- data.frame(metric.name, metric.val)
options(digits = 4)
D

# Linear Regression Line
yhat <- lm(Y ~ X1)
abline(yhat, col="steelblue", lwd = 3)

residual <- resid(yhat)
plot(X1, residual, xlab = "Deaths", ylab = "residuals for model",
     main="Residual plot for our linear model")
abline(0, 0, col ="red", lwd=3)

alpha<-0.05
tcrits<-qt(alpha/2, df=n-2, lower.tail = F)
beta1<-0
epsS <- tcrits*SE.beta1hat
tstats <- (beta1hat - beta1)/SE.beta1hat
CIL <- beta1hat - epsS
CIU <- beta1hat + epsS
pvalS <- 2*pt(abs(tstats), df=n-2, lower.tail = F)
pvalS


#Part II.B – Simple Quadratic Regression
M2 <- nemolm2(Y, cbind(X1, X1^2))

#Standardized residual plot
plot(X1, M2$residual, xlab = "Crashes", ylab = "residuals for Simple Quadratic model",
     main="Residual plot for Simple Quadratic model")
abline(0, 0, col ="red", lwd=3)
M2

#Part III – Multiple Linear Regression
M3 <- nemolm2(Y, cbind(X1, X2, X3))
M3

#ANOVA table
metric_name_A <-c("SST", "MST", "SSM", "MSM", "SSE", "MSE", "Fstat", "p-value")
metric_val_A <-c(M3$sst, M3$mst, M3$ssm, M3$msm, M3$sse, M3$mse, M3$fstat, M3$pval)
Summary_A <- data.frame(metric_name_A, metric_val_A)
Summary_A

#Variance inflation factors calculated for each variable with barplot
# Y regressed on X1, X2, and X3
MYvX1c <- nemolm2(Y, cbind(X2, X3))
MX1vX1c <- nemolm2(X1, cbind(X2, X3))

MYvX2c <- nemolm2(Y, cbind(X1, X3))
MX2vX2c <- nemolm2(X2, cbind(X1, X3))

MYvX3c <- nemolm2(Y, cbind(X1, X2))
MX3vX3c <- nemolm2(X3, cbind(X1, X2))

vif1 <- 1/(1-MYvX1c$r2)
vif2 <- 1/(1-MYvX2c$r2)
vif3 <- 1/(1-MYvX3c$r2)

vif <- c(vif1, vif2, vif3)
barplot(vif, horiz=T, main="Variance Inflation Factors", 
        names.arg = c('Crashes', 'Miles traveled (millions)', 'Motor vehicles'), 
        xlim=c(0,1425))

#new fits
M4 <- nemolm2(Y, cbind(X1, X2, X3, X2*X3))
M4

#Added variable plots for each variable
plot(MYvX1c$sres, MX1vX1c$sres, 
     main="Added Variable Plot fot X1", 
     xlab = "S.Residuals for Y~X1c", 
     ylab = "S.Residuals for X1~X1c")
abline(0,0, lwd=2)
abline(mean(MX1vX1c$sres)-cor(MYvX1c$sres, 
                              MX1vX1c$sres)*sd(MX1vX1c$sres)/sd(MYvX1c$sres)*mean(MYvX1c$sres),
       cor(MYvX1c$sres, MX1vX1c$sres)*sd(MX1vX1c$sres)/sd(MYvX1c$sres))

plot(MYvX2c$sres, MX2vX2c$sres, 
     main="Added Variable Plot fot X2", 
     xlab = "S.Residuals for Y~X2c",
     ylab = "S.Residuals for X2~X2c")
abline(0,0, lwd=2)
abline(mean(MX2vX2c$sres)-cor(MYvX2c$sres, 
                              MX2vX2c$sres)*sd(MX2vX2c$sres)/sd(MYvX2c$sres)*mean(MYvX2c$sres), 
       cor(MYvX2c$sres, MX2vX2c$sres)*sd(MX2vX2c$sres)/sd(MYvX2c$sres))

plot(MYvX3c$sres, MX3vX3c$sres, 
     main="Added Variable Plot fot X3", 
     xlab = "S.Residuals for Y~X3c", 
     ylab = "S.Residuals for X3~X3c")
abline(0,0, lwd=2)
abline(mean(MX3vX3c$sres)-cor(MYvX3c$sres, 
                              MX3vX3c$sres)*sd(MX3vX3c$sres)/sd(MYvX3c$sres)*mean(MYvX3c$sres), 
       cor(MYvX3c$sres, MX3vX3c$sres)*sd(MX3vX3c$sres)/sd(MYvX3c$sres))

#Standardized residual plot with title and axis labels
plot(X1, M2$std.residual, xlab = "Year", ylab = "residuals for model",
     main="Residual plot for Simple Quadratic model")
abline(0, 0, col ="red", lwd=3)

#Construction of the correlation matrix between Y and all three variables
cor(S)

#Part IV – Time Series Fundamentals
plot(S$Year, S$Deaths, type='o')

#Standardize
XT1 <- S$Year
YT1 <- S$Deaths

sxT1  <- sd(XT1)
xT1 <- mean(XT1)
XT <- (XT1-xT1)/sxT1

syT1  <- sd(YT1)
yT1 <- mean(YT1)
YT <- (YT1-yT1)/syT1

MT<-nemolm2(YT, cbind(XT, XT^2, XT^3, XT^4))
MT
plot(XT, YT, type='o')
lines(XT, MT$predicted, lwd=3, col='green')
