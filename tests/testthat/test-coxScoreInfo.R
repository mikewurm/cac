survivalInst <- require(survival)

## Test cacInfoApprox functions ###############################################
set.seed(100)
n <- 100
beta <- -2:2
p <- length(beta)
x <- matrix(rnorm(n*p), nrow=n)
hr <- exp(c(x %*% beta))
y <- rexp(n, rate=hr)
d <- rbinom(n, 1, .8)

test_that("coxScoreInfo matches coxph", {
    if (!survivalInst) skip("survival not installed.")

    # Fit coxph with arbitrary weights and offset
    wts <- runif(n)
    off <- rnorm(n)
    fit <- coxph(Surv(y, d) ~ x + offset(off), weights=wts)
    si <- cac:::coxProfileScoreInfo(beta=coef(fit), x, y, d, wt=wts, offset=off)

    expect_equal(si$profileScore, rep(0, p))
    expect_equal(solve(si$profileInfo), fit$var)
})

# ## Verify that weights affect the info
# fit2 <- coxph(Surv(y, d) ~ x + offset(off))  # without weights
# all.equal(solve(si$info), fit2$var)  # should be different
#
# ## Verify that offset affects the info
# fit3 <- coxph(Surv(y, d) ~ x, weights=wts)  # without offset
# all.equal(solve(si$info), fit3$var)  # should be different
