library(survival)

## Test cacInfoApprox functions ###############################################
set.seed(100)
n <- 100
beta <- cbind(-2:2, 2:-2, rep(0, 5))
pi <- c(.25, .25, .5)
p <- nrow(beta)
k <- ncol(beta)
u <- sample(1:k, n, replace=TRUE, prob=pi)
x <- matrix(rnorm(n*p), nrow=n)
hr <- exp((x %*% beta)[cbind(1:n, u)])
y <- rexp(n, rate=hr)
d <- rbinom(n, 1, .8)

## True baseline hazard functions are identical
uStart <- matrix(0, nrow=n, ncol=k)
uStart[cbind(1:n, u)] <- 1
m1 <- cac(x, y, d, k=3, hazRestrict="identical", uStart=list(uStart))
m2 <- cac(x, y, d, k=3, hazRestrict="proportional", uStart=list(uStart))
m3 <- cac(x, y, d, k=3, hazRestrict="none", uStart=list(uStart))

## Test coxScoreInfo
wt <- runif(n, 0, 5)
coxmod <- coxph(survival::Surv(y, d) ~ x, weights=wt, ties="breslow")
si <- coxScoreInfo(coxmod$coefficients, x, y, d, wt)
max(abs(si$score)) < 1e-7
all.equal(solve(si$info), coxmod$var)

## NEW: compare with "terms" function
si2 <- coxScoreInfo2(coxmod$coefficients, x, y, d, wt)
all.equal(si$score, si2$score)
all.equal(si$info, si2$info)
microbenchmark::microbenchmark(
    coxScoreInfo(coxmod$coefficients, x, y, d, wt),
    coxScoreInfo2(coxmod$coefficients, x, y, d, wt),
    times=100
)

## Test cacSI.huProfiled
a1 <- cacSI.huProfiled(m1)
a2 <- cacSI.huProfiled(m2)
a3 <- cacSI.huProfiled(m3)
max(abs(c(a1$piScore, a1$betaScore))) < 1e-7
max(abs(c(a2$piScore, a2$betaScore))) < 1e-7
max(abs(c(a3$piScore, a3$betaScore))) < 1e-7

## Test cacMSE with (pi, beta) at the MLE: should be the same as cacSI.huProfiled
b1 <- cacMSI(m1$betaHat, m1$piHat, m1)
b2 <- cacMSI(m2$betaHat, m2$piHat, m2)
b3 <- cacMSI(m3$betaHat, m3$piHat, m3)
all.equal(c(b1$piM, b1$betaM), c(m1$piHat[-k], m1$betaHat))
all.equal(c(b2$piM, b2$betaM), c(m2$piHat[-k], m2$betaHat))
all.equal(c(b3$piM, b3$betaM), c(m3$piHat[-k], m3$betaHat))
max(abs(c(b1$piScore, b1$betaScore))) < 1e-7
max(abs(c(b2$piScore, b2$betaScore))) < 1e-7
max(abs(c(b3$piScore, b3$betaScore))) < 1e-7
all.equal(list(b1$piPECDInfo, b1$betaPECDInfo), list(a1$piPECDInfo, a1$betaPECDInfo))
all.equal(list(b2$piPECDInfo, b2$betaPECDInfo), list(a2$piPECDInfo, a2$betaPECDInfo))
all.equal(list(b3$piPECDInfo, b3$betaPECDInfo), list(a3$piPECDInfo, a3$betaPECDInfo))

## Test cacMSI with (pi, beta) at true values, rather than MLE (just make sure it runs)
c1 <- cacMSI(beta, pi, m1)
c2 <- cacMSI(beta, pi, m2)
c3 <- cacMSI(beta, pi, m3)

## Test cacInfoApprox (make sure it runs, see if pi variance is close to louis version)
d1 <- cacInfoApprox(m1)
d2 <- cacInfoApprox(m2)
d3 <- cacInfoApprox(m3)

## NEW: compare with louisProfiledInfo
louisInfo(m1)

e1 <- louisInfoPi(m1)
e2 <- louisInfoPi(m2)
e3 <- louisInfoPi(m3)

d1$varApproxM[1:(k-1), 1:(k-1)]
solve(e1)
d2$varApproxM[1:(k-1), 1:(k-1)]
solve(e2)
d3$varApproxM[1:(k-1), 1:(k-1)]
solve(e3)
