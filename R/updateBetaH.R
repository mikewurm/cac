## Interesting note about functions that accept formula:
## y~as.matrix(x) works faster
## y~., data=as.data.frame(x) works but slower
## y~as.data.frame(x) does not work
## y~., data=as.matrix(x) does not work

## offsetMat is only used if !is.null(fixedBeta)
## fixedInt is only used if hazRestrict="proportional" and !is.null(fixedBeta)
updateBetaH <- function(hazRestrict, coxX, coxY, coxInt, uniqY, yCt, uHat,
                        fixedBeta, offsetMat)
{
    uHat <- uHat + 1e-100  # to prevent coxph "Error during wrapup: NA/NaN/Inf in foreign function call (arg 6)"
    # posCls <- which(apply(uHat, 2, function(x) any(x > 1e-100)))
    n <- nrow(uHat)
    k <- ncol(uHat)
    p <- ncol(coxX) / k^(hazRestrict!="none")
    betaHat <- matrix(0, nrow=p, ncol=k)
    HHatUniq <- matrix(0, nrow=length(uniqY), ncol=k)

    if (hazRestrict=="none")
    {
        # for (i in posCls)
        for (i in 1:k)
        {
            if (is.null(fixedBeta))
            {
                # mod <- suppressWarnings(survival::coxph(coxY ~ coxX,
                #     weights=uHat[,i], subset=uHat[,i]>0, ties="breslow"))
                mod <- suppressWarnings(survival::coxph(coxY ~ coxX, weights=uHat[,i], ties="breslow"))
                coeff <- coef(mod)
                coeff[is.na(coeff)] <- 0  # in case of singular X
                betaHat[, i] <- coeff
                blHaz <- suppressWarnings(survival::basehaz(mod, centered=FALSE))
                HHatUniq[, i] <- if (nrow(blHaz) == length(uniqY)) {
                    blHaz$hazard
                } else {
                    getH(uniqY, blHaz)
                }
            } else
            {
                betaHat <- fixedBeta
                # mod <- suppressWarnings(survival::coxph(coxY ~ offset(offsetMat[, i]),
                #     weights=uHat[,i], subset=uHat[,i]>0, ties="breslow"))
                mod <- suppressWarnings(survival::coxph(coxY ~ offset(offsetMat[, i]), weights=uHat[,i], ties="breslow"))
                blHaz <- suppressWarnings(survival::basehaz(mod, centered=TRUE))
                    # centered=FALSE does not work with offset-only model
                HHatUniq[, i] <- if (nrow(blHaz) == length(uniqY)) {
                    blHaz$hazard
                } else {
                    getH(uniqY, blHaz)
                }
                HHatUniq[, i] <- HHatUniq[, i] / exp(mean(offsetMat[uHat[, i]>0, i]))
                    # centered=FALSE does not work with offset-only model
            }
        }
    } else if (hazRestrict=="proportional")
    {
        # keepRow <- do.call(c, lapply(posCls, function(x) 1:n + n*(x-1)))
        # coxY <- coxY[keepRow, ]
        # keepX <- do.call(c, lapply(posCls, function(x) 1:p + p*(x-1)))
        # coxX <- coxX[keepRow, keepX]
        # uHat <- as.vector(uHat[, posCls])
        # coxInt <- coxInt[keepRow, posCls]
        uHat <- c(uHat)
        if (is.null(fixedBeta)) {
            # mod <- suppressWarnings(survival::coxph(coxY ~ cbind(coxInt, coxX),
            #     weights=uHat, subset=uHat>0, ties="breslow"))
            mod <- suppressWarnings(survival::coxph(coxY ~ cbind(coxInt, coxX), weights=uHat, ties="breslow"))
            coeff <- coef(mod)
            coeff[is.na(coeff)] <- 0  # in case of singular X
            # betaHat[, posCls] <- coeff[-(1:length(posCls))]
            betaHat[] <- coeff[-(1:k)]
            # basehazFact <- rep(1, k)
            # basehazFact[posCls] <- exp(coeff[1:length(posCls)])
            # basehazFact <- exp(coeff[1:length(posCls)])
            basehazFact <- exp(coeff[1:k])
            blHaz <- suppressWarnings(survival::basehaz(mod, centered=FALSE))
            HHatUniq <- if (nrow(blHaz) == length(uniqY)) {
                blHaz$hazard
            } else {
                getH(uniqY, blHaz)
            }
            HHatUniq <- HHatUniq * matrix(rep(basehazFact, each=length(uniqY)), nrow=length(uniqY))
                # converts HHatUniq from vector to matrix
        } else {
            intHat <- fixedBeta[1, ]
            betaHat <- fixedBeta[-1, , drop=FALSE]
            # basehazFact <- rep(1, k)
            # basehazFact[posCls] <- exp(intHat[posCls])
            basehazFact <- exp(intHat)
            # offsetVec <- as.vector(offsetMat[, posCls])
            offsetVec <- c(offsetMat)
            # mod <- suppressWarnings(survival::coxph(coxY ~ offset(offsetVec),
            #     weights=uHat, subset=uHat>0, ties="breslow"))
            mod <- suppressWarnings(survival::coxph(coxY ~ offset(offsetVec), weights=uHat, ties="breslow"))
            blHaz <- suppressWarnings(survival::basehaz(mod, centered=TRUE))
                # centered=FALSE does not work with offset-only model
            HHatUniq <- if (nrow(blHaz) == length(uniqY)) {
                blHaz$hazard
            } else {
                getH(uniqY, blHaz)
            }
            # HHatUniq <- HHatUniq / exp(mean(offsetVec[uHat>0])) * matrix(basehazFact, nrow=length(uniqY), ncol=k, byrow=TRUE)
            HHatUniq <- HHatUniq / exp(mean(offsetVec)) * matrix(basehazFact, nrow=length(uniqY), ncol=k, byrow=TRUE)
                # centered=FALSE does not work with offset-only model
        }
    } else if (hazRestrict=="identical")
    {
        # keepRow <- do.call(c, lapply(posCls, function(x) 1:n + n*(x-1)))
        # coxY <- coxY[keepRow, ]
        # keepX <- do.call(c, lapply(posCls, function(x) 1:p + p*(x-1)))
        # coxX <- coxX[keepRow, keepX]
        # uHat <- as.vector(uHat[, posCls])
        uHat <- c(uHat)
        if (is.null(fixedBeta)) {
            # mod <- suppressWarnings(survival::coxph(coxY ~ coxX,
            #     weights=uHat, subset=uHat>0, ties="breslow"))
            mod <- suppressWarnings(survival::coxph(coxY ~ coxX, weights=uHat, ties="breslow"))
            coeff <- coef(mod)
            coeff[is.na(coeff)] <- 0  # in case of singular X
            # betaHat[, posCls] <- coeff
            betaHat[] <- coeff
            blHaz <- suppressWarnings(survival::basehaz(mod, centered=FALSE))
            HHatUniq <- if (nrow(blHaz) == length(uniqY)) {
                blHaz$hazard
            } else {
                getH(uniqY, blHaz)
            }
            HHatUniq <- matrix(rep(HHatUniq, k), ncol=k)
        } else {
            betaHat <- fixedBeta
            # offsetVec <- as.vector(offsetMat[, posCls])
            offsetVec <- c(offsetMat)
            # mod <- suppressWarnings(survival::coxph(coxY ~ offset(offsetVec),
            #     weights=uHat, subset=uHat>0, ties="breslow"))
            mod <- suppressWarnings(survival::coxph(coxY ~ offset(offsetVec), weights=uHat, ties="breslow"))
            blHaz <- suppressWarnings(survival::basehaz(mod, centered=TRUE))
                # centered=FALSE does not work with offset-only model
            HHatUniq <- if (nrow(blHaz) == length(uniqY)) {
                blHaz$hazard
            } else {
                getH(uniqY, blHaz)
            }
            # HHatUniq <- matrix(rep(HHatUniq / exp(mean(offsetVec[uHat>0])), k), ncol=k)
            HHatUniq <- matrix(rep(HHatUniq / exp(mean(offsetVec)), k), ncol=k)
                # centered=FALSE does not work with offset-only model
        }
    }

    hHatUniq <- apply(HHatUniq, 2, function(x) x - c(0, head(x, -1)))
    hHat <- apply(hHatUniq, 2, rep, yCt)
    HHat <- apply(HHatUniq, 2, rep, yCt)
    list(betaHat=betaHat, hHat=hHat, HHat=HHat)
}

getH <- function(uniqY, blHaz)
{
    timeVals <- c(blHaz$time, Inf)
    Hvals <- c(0, blHaz$hazard)
    index <- 1
    timeVal <- timeVals[index]
    Hval <- Hvals[index]
    sapply(uniqY, function(x) {
        while (x >= timeVal) {
            index <<- index+1
            Hval <<- Hvals[index]
            timeVal <<- timeVals[index]
        }
        Hval
    })
}

# getBasehaz <- function(coxMod)
# {
#     bh <- suppressWarnings(survival::basehaz(coxMod, centered=FALSE))
#     basehaz <- with(bh, data.frame(
#         time=time,
#         hHat=c(hazard[1], hazard[-1] - hazard[-nrow(bh)]),
#         HHat=hazard))
#     basehaz
# }

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
# betaH1 <- updateBetaH("none", coxX, coxY, coxInt=NULL, uniqY, yCt, uHat, piHat, fixedBeta=NULL)
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
# betaH3$hHatUniq[,4:5] / betaH3$hHatUniq[,3]
# betaH4$hHatUniq[,4:5] / betaH4$hHatUniq[,3]
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

