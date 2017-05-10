## Interesting note about functions that accept formula:
## y~as.matrix(x) works faster
## y~., data=as.data.frame(x) works but slower
## y~as.data.frame(x) does not work
## y~., data=as.matrix(x) does not work
## penalized package is not an option because it does not allow regression weights

## xBeta is only used if !is.null(fixedBeta)
updateBetaH <- function(hazRestrict, coxX, coxY, coxInt, yCt, uHat, piHat, fixedBeta, xBeta)
{
    n <- nrow(uHat)
    k <- ncol(uHat)
    p <- ncol(coxX) / k^(hazRestrict!="none")
    uHat[uHat==0] <-1e-323 # coxph weights must be >0
    # uHat[uHat==0] <- min(uHat[uHat>0]) * 1e-4 # coxph weights must be >0
    posCls <- which(piHat != 0)
    betaHat <- matrix(0, nrow=p, ncol=k)
    hHatOrd <- matrix(0, nrow=length(unique(coxY[,1])), ncol=k)

    if (hazRestrict=="none")
    {
        for (i in posCls)
        {
            if (is.null(fixedBeta)) {
                mod <- suppressWarnings(survival::coxph(coxY ~ coxX,
                    weights=uHat[,i], ties="breslow"))
                coeff <- coef(mod)
                coeff[is.na(coeff)] <- 0  # in case of singular X
                betaHat[, i] <- coeff
                blHaz <- getBasehaz(mod)
                hHatOrd[, i] <- blHaz$hHat
            } else {
                betaHat <- fixedBeta
                mod <- suppressWarnings(survival::coxph(coxY ~ offset(xBeta[, i]),
                    weights=uHat[,i], ties="breslow"))
                blHaz <- getBasehaz(mod)
                # centered=FALSE does not work for basehaz with offset-only model
                hHatOrd[, i] <- blHaz$hHat / exp(mean(xBeta[, i]))
            }
        }
    } else if (hazRestrict=="proportional")
    {
        keepRow <- do.call(c, lapply(posCls, function(x) 1:n + n*(x-1)))
        coxY <- coxY[keepRow, ]
        keepX <- do.call(c, lapply(posCls, function(x) 1:p + p*(x-1)))
        coxX <- coxX[keepRow, keepX]
        uHat <- as.vector(uHat[, posCls])
        coxInt <- coxInt[keepRow, posCls]
        if (is.null(fixedBeta)) {
            mod <- suppressWarnings(survival::coxph(coxY ~ cbind(coxInt, coxX),
                weights=uHat, ties="breslow"))
            coeff <- coef(mod)
            coeff[is.na(coeff)] <- 0  # in case of singular X
            betaHat[, posCls] <- coeff[-(1:length(posCls))]
            basehazFact <- rep(0, k)
            basehazFact[posCls] <- exp(coeff[1:length(posCls)])
            blHaz <- getBasehaz(mod)
            hHatOrd <- blHaz$hHat * matrix(rep(basehazFact, each=nrow(blHaz)), nrow=nrow(blHaz))
        } else {
            betaHat <- fixedBeta
            xBeta <- as.vector(xBeta[, posCls])
            mod <- suppressWarnings(survival::coxph(coxY ~ offset(xBeta) + coxInt,
                weights=uHat, ties="breslow"))
            coeff <- coef(mod)
            coeff[is.na(coeff)] <- 0  # in case of singular X
            basehazFact <- rep(0, k)
            basehazFact[posCls] <- exp(coeff)
            blHaz <- getBasehaz(mod)
            hHatOrd <- blHaz$hHat * matrix(rep(basehazFact, each=nrow(blHaz)), nrow=nrow(blHaz))
        }
    } else if (hazRestrict=="identical")
    {
        keepRow <- do.call(c, lapply(posCls, function(x) 1:n + n*(x-1)))
        coxY <- coxY[keepRow, ]
        keepX <- do.call(c, lapply(posCls, function(x) 1:p + p*(x-1)))
        coxX <- coxX[keepRow, keepX]
        uHat <- as.vector(uHat[, posCls])
        if (is.null(fixedBeta)) {
            mod <- suppressWarnings(survival::coxph(coxY ~ coxX,
                weights=uHat, ties="breslow"))
            coeff <- coef(mod)
            coeff[is.na(coeff)] <- 0  # in case of singular X
            betaHat[, posCls] <- coeff
            blHaz <- getBasehaz(mod)
            hHatOrd <- matrix(rep(blHaz$hHat, k), ncol=k)
        } else {
            betaHat <- fixedBeta
            xBeta <- as.vector(xBeta[, posCls])
            mod <- suppressWarnings(survival::coxph(coxY ~ offset(xBeta),
                                                    weights=uHat, ties="breslow"))
            blHaz <- getBasehaz(mod)
            hHatOrd <- matrix(rep(blHaz$hHat / exp(mean(xBeta)), k), ncol=k)
        }
    }

    HHatOrd <- apply(hHatOrd, 2, cumsum)
    hHat <- apply(hHatOrd, 2, rep, yCt)
    HHat <- apply(HHatOrd, 2, rep, yCt)
    hHatOrd <- data.frame(hHatOrd)
    HHatOrd <- data.frame(HHatOrd)
    names(hHatOrd) <- paste0("h", 1:k)
    names(HHatOrd) <- paste0("H", 1:k)
    hHatOrd <- data.frame(time=blHaz$time, hHatOrd)
    HHatOrd <- data.frame(time=blHaz$time, HHatOrd)
    list(betaHat=betaHat, hHat=hHat, HHat=HHat, hHatOrd=hHatOrd, HHatOrd=HHatOrd)
}

getBasehaz <- function(coxMod)
{
    bh <- suppressWarnings(survival::basehaz(coxMod, centered=FALSE))
    basehaz <- with(bh, data.frame(
        time=time,
        hHat=c(hazard[1], hazard[-1] - hazard[-nrow(bh)]),
        HHat=hazard))
    basehaz
}

## Not needed?
# tryCatch(
# {
#     xx <- coxX[posWtID, ]
#     yy <- coxY[posWtID]
#     ww <- uHat[posWtID, k]
#     mod <- survival::coxph(yy ~ xx, weights=ww, ties="breslow")
#     coef(mod) # coef does not appear to be a method from the survival package
# }, error=function(cond)
# {
#     warning("coxph error.")
#     rep(NA, p)
# })

## Test - includes tied times and zero weights (uHat=0)
# set.seed(123)
# n <- 100
# p <- 5
# piHat <- c(0, c(1, 1, 1) / 3)
# pCens <- .5
#
# k <- length(piHat)
# uHat <- cbind(0, 1, matrix(rpois(n*(k-2), 1), nrow=n))
# uHat <- uHat / rowSums(uHat)
#
# coxX <- matrix(rnorm(n*p), nrow=n)
# y <- rep(rexp(n/2), 2)
# d <- rbinom(n, 1, pCens)
# coxY <- survival::Surv(y, d)
# uniqY <- unique(y)
# yCt <- sapply(uniqY, function(x) sum(y==x))
# fixedBeta <- matrix(rnorm(p*k), nrow=p)
# xBeta <- coxX%*%fixedBeta
#
# makeBlockDiag <- function(x, p, k)
# {
#     bdRow <- function(i) {
#         leftBlock <- c(rep(list(0), p*(i-1)))
#         midBlock <- list(x)
#         rightBlock <- rep(list(0), p*(k-i))
#         blockList <- c(leftBlock, midBlock, rightBlock)
#         do.call(cbind, blockList)
#     }
#     blockX <- do.call(rbind, lapply(1:k, bdRow))
#     blockX
# }
# coxY2 <- survival::Surv(rep(y, k), rep(d, k))
# coxX2 <- makeBlockDiag(coxX, p, k)
# coxInt <- sapply(1:k, function(i) c(rep(0, n*(i-1)), rep(1, n), rep(0, n*(k-i))))
#
# betaH1 <- updateBetaH("none", coxX, coxY, coxInt=NULL, yCt, uHat, piHat, fixedBeta=NULL)
# betaH1
# betaH2 <- updateBetaH("none", coxX, coxY, coxInt=NULL, yCt, uHat, piHat, fixedBeta, xBeta)
# betaH2
# betaH3 <- updateBetaH("proportional", coxX2, coxY2, coxInt, yCt, uHat, piHat, fixedBeta=NULL)
# betaH3
# betaH4 <- updateBetaH("proportional", coxX2, coxY2, coxInt, yCt, uHat, piHat, fixedBeta, xBeta)
# betaH4
# betaH5 <- updateBetaH("identical", coxX2, coxY2, coxInt=NULL, yCt, uHat, piHat, fixedBeta=NULL)
# betaH5
# betaH6 <- updateBetaH("identical", coxX2, coxY2, coxInt=NULL, yCt, uHat, piHat, fixedBeta, xBeta)
# betaH6
# betaHList <- list(betaH1, betaH2, betaH3, betaH4, betaH5, betaH6)
# sapply(betaHList, function(x) x$betaHat)
# betaH3$hHatOrd[,4:5] / betaH3$hHatOrd[,3]
# betaH4$hHatOrd[,4:5] / betaH4$hHatOrd[,3]
#
# ## test: coxph cannot fit a model with all weights at 1e-323
# t1 <- survival::coxph(coxY ~ coxX, weights=rep(1e-323, n))
#
# ## test: zero weights are fine with offset
# t2 <- survival::coxph(coxY ~ offset(rnorm(n)), weights=rep(1e-323, n))
#
# ## test: zero weights are fine with offset
# t3 <- survival::coxph(coxY ~ offset(rnorm(n)), weights=runif(n))
#
# ## test: weights do affect baseline hazard estimation (as expected)
# t4 <- survival::coxph(coxY ~ offset(rnorm(n)), weights=c(rep(c(1e-323, 1), each=n/2)))
# do.call(cbind, lapply(list(t2, t3, t4), function(x) basehaz(x)$hazard))
#
# ## test
# yy <- survival::Surv(rep(y, k), rep(d, k))
# xx <- sapply(2:k, function(i) c(rep(0, n*(i-1)), rep(1, n), rep(0, n*(k-i))))
# t5 <- survival::coxph(yy~xx)
#
# xx2 <- sapply(1:k, function(i) c(rep(0, n*(i-1)), rep(1, n), rep(0, n*(k-i))))
# t6 <- survival::coxph(yy~xx2)
# cbind(coef(t5), coef(t6))
#
# mod <- suppressWarnings(survival::coxph(coxY ~ cbind(coxInt, coxX),
#                                         weights=uHat, ties="breslow"))
