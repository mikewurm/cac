# Approximates the observed information of the observed data log-likelihood
# Calculates both the NDM and NDS approximations using Richardson extrapolation

# Note that for hazRestrict=="proportional", the intercepts are included as the first k-1 terms of beta

#' CAC profile information
#'
#' Calculates standard errors for pi and beta (and intercepts under the proportional
#' baseline hazard model). The observed information of the semiparametric likelihood
#' is approximated by two numerical differentiation techniques: NDM and NDS.
#' First order Richardson extrapolation is used for both techniques.
#' The NDM approach differentiates the EM algorithm operator, and the NDS approach
#' numerically differentiates the profile score function. These methods should return
#' similar approximations.
#'
#' @param cacFit An S3 object of class 'cacFit' returned from the cac function.
#' @param eta Determines the step size around the MLE for numerical differentiation.
#'
#' @return A list containing the following:
#' \describe{
#'   \item{sePi}{Vector of standard errors corresponding to the mixture probability estimates.}
#'   \item{seBeta}{Matrix of standard errors corresponding to the coefficient estimates.}
#'   \item{seInt}{Vector of standard errors corresponding to the intercept estimates, which
#'   represent the log hazard ratio between the baseline hazard functions. Only applicable
#'   when \code{hazRestrict = "proportional"}.}
#'   \item{var}{Variance matrix estimate for pi, intercepts (if applicable), and beta,
#'   in that order.}
#'   \item{profileInfo}{The approximate observed information matrix. Has fewer components
#'   than the variance matrix because baseline classes for pi and intercept (if applicable)
#'   estimators are not included, because these would make the matrix non-invertible.}
#'   \item{profileInfoRaw}{The approximate observed information matrix, before adjusting
#'   for symmetry. \code{profileInfo = (profileInfo + t(profileInfo)) / 2}.}
#' }
#' Each item has NDM and NDS versions, corresponding to the two approximation methods.
#'
#' @export
cacProfileInfo <- function(cacFit, eta=1e-4)
{
    betaHat <- cacFit$betaHat
    piHat <- cacFit$piHat
    p <- nrow(betaHat)
    k <- ncol(betaHat)
    hazRestrict <- cacFit$arglist$hazRestrict

    nPar <- k-1 + k*p  # number of parameters in profile likelihood
    if (hazRestrict == "proportional") {
        nPar <- nPar + k-1  # include intercepts
        hHat <- cacFit$hHat
        intHat <- colMeans(log(hHat / hHat[, k]), na.rm=TRUE)
    } else {
        intHat <- NULL
    }
    dm <- dscore <- matrix(NA, nrow=nPar, ncol=nPar)
    mList <- scoreList <- rep(list(matrix(NA, nrow=nPar, ncol=nPar)), 4)
    h <- rep(NA, nPar)
    reFact <- c(-2, -1, 1, 2)  # Richardson extrapolation factors

    for (i in 1:nPar) {
        delPi <- rep(0, k)
        delBeta <- rep(0, nPar-(k-1))
        if (i < k) {
            delPi[i] <- h[i] <- eta
            delPi[k] <- -h[i]
        } else {
            delBeta[i-(k-1)] <- h[i] <- eta * max(abs(c(betaHat, intHat[-k])[i-(k-1)]), 1)
        }

        for (q in 1:4) {
            newpi <- piHat + reFact[q] * delPi
            if (hazRestrict == "proportional") {
                newbeta <- betaHat + reFact[q] * head(delBeta, -(k-1))
                newint <- intHat + c(reFact[q] * tail(delBeta, k-1), 0)
                    # last intercept is fixed at zero and is not a parameter
                newbeta <- rbind(newint, newbeta)
            } else {
                newbeta <- betaHat + reFact[q] * delBeta
            }
            msi <- cacProfileMSI(pi=newpi, beta=newbeta, cacFit=cacFit)
            mList[[q]][, i] <- with(msi, c(mPi, mBeta))
                # mPi has length k-1, mBeta is a vector of length k*p (or k-1+k*p if hazRestrict="proportional)
            scoreList[[q]][, i] <- with(msi, c(profileScorePi, profileScoreBeta))
        }

    }

    dm <- (mList[[1]] - 8*mList[[2]] + 8*mList[[3]] - mList[[4]]) / (12*rep(h, each=nPar))
    dscore <- (scoreList[[1]] - 8*scoreList[[2]] + 8*scoreList[[3]] - scoreList[[4]]) / (12*rep(h, each=nPar))

    # Get expected complete data profile information for NDM approximation
    si <- cacProfileScoreInfoECD(cacFit)

    # NDM and NDS profile info approximations
    profileInfoECDPi <- si$profileInfoECDPi
    profileInfoECDBeta <- si$profileInfoECDBeta
    profileInfoECD <- matrix(0, nrow=nPar, ncol=nPar)
    profileInfoECD[1:(k-1), 1:(k-1)] <- profileInfoECDPi
    profileInfoECD[k:nPar, k:nPar] <- profileInfoECDBeta
    profileInfoRaw.NDM <- profileInfoECD %*% (diag(nPar) - dm)
    profileInfoRaw.NDS <- -dscore

    ## Force the approximate information to be symmetric
    profileInfo.NDM <- (profileInfoRaw.NDM + t(profileInfoRaw.NDM)) / 2
    profileInfo.NDS <- (profileInfoRaw.NDS + t(profileInfoRaw.NDS)) / 2

    ## Get standard errors
    if (hazRestrict == "proportional") {
        lincom1 <- diag(nPar)
        lincom <- rbind(lincom1[1:(k-1), ], rep(c(-1, 0), c(k-1, nPar-(k-1))), lincom1[k:nPar, ], rep(0, nPar))
        var.NDM <- lincom %*% solve(profileInfo.NDM) %*% t(lincom)
        var.NDS <- lincom %*% solve(profileInfo.NDS) %*% t(lincom)
        seAll.NDM <- sqrt(diag(var.NDM))
        seAll.NDS <- sqrt(diag(var.NDS))
        sePi.NDM <- seAll.NDM[1:k]
        sePi.NDS <- seAll.NDS[1:k]
        seBeta.NDM <- matrix(seAll.NDM[seq(k+1, length.out=k*p)], ncol=k)
        seBeta.NDS <- matrix(seAll.NDS[seq(k+1, length.out=k*p)], ncol=k)
        seInt.NDM <- tail(seAll.NDM, k)
        seInt.NDS <- tail(seAll.NDS, k)
    } else {
        lincom1 <- diag(nPar)
        lincom <- rbind(lincom1[1:(k-1), ], rep(c(-1, 0), c(k-1, nPar-(k-1))), lincom1[k:nPar, ])
        var.NDM <- lincom %*% solve(profileInfo.NDM) %*% t(lincom)
        var.NDS <- lincom %*% solve(profileInfo.NDS) %*% t(lincom)
        seAll.NDM <- sqrt(diag(var.NDM))
        seAll.NDS <- sqrt(diag(var.NDS))
        sePi.NDM <- seAll.NDM[1:k]
        sePi.NDS <- seAll.NDS[1:k]
        seBeta.NDM <- matrix(seAll.NDM[seq(k+1, length.out=k*p)], ncol=k)
        seBeta.NDS <- matrix(seAll.NDS[seq(k+1, length.out=k*p)], ncol=k)
        seInt.NDM <- NA
        seInt.NDS <- NA
    }

    list(sePi.NDM=sePi.NDM, sePi.NDS=sePi.NDS,
         seBeta.NDM=seBeta.NDM, seBeta.NDS=seBeta.NDS,
         seInt.NDM=seInt.NDM, seInt.NDS=seInt.NDS,
         var.NDM=var.NDM, var.NDS=var.NDS,
         profileInfo.NDM=profileInfo.NDM, profileInfo.NDS=profileInfo.NDS,
         profileInfoRaw.NDM=profileInfoRaw.NDM, profileInfoRaw.NDS=profileInfoRaw.NDS)
}

# Computes EM operator M(.), profile score, and expected profile observed info.
# beta New beta (either beta or pi should be a small step away from the cacFit estimate).
# pi New pi (full vector, no baseline class) (either beta or pi should be a small
#    step away from the cacFit estimate).
# cacFit S3 object of class 'cacFit'. The data set (x, y, d) is extracted from here,
#        as well as starting values for uHat. Random starts are not needed if the new
#        (beta, pi) are close to the cacFit estimate, as this function is intended.
# Returns mBeta as a vector, not a matrix. If hazRestrict="parallel", then the k-1 (not k) intercepts are first.
## This function is only useful for calculating the score near a local CAC solution.
## hHat and HHat must be provided. (Otherwise random starts would be necessary to
## find uHat and hHat/HHat).
cacProfileMSI <- function(pi, beta, cacFit)
{
    arglist <- cacFit$arglist
    n <- nrow(arglist$x)
    p <- ncol(arglist$x)
    k <- ncol(cacFit$betaHat)
    x <- arglist$x
    y <- arglist$y
    d <- arglist$d
    maxiter <- arglist$maxiter
    eps <- arglist$eps
    hazRestrict <- arglist$hazRestrict
    betaHat <- cacFit$betaHat
    uHat <- cacFit$uHat

    ## Profile over baseline hazard (update hHat, HHat and uHat with beta and pi fixed)
    cacFitProfiled <- cac(x=x, y=y, d=d, k=k, hazRestrict=hazRestrict,
                          maxiter=maxiter, eps=eps, uStart=list(uHat),
                          fixedBeta=beta, fixedPi=pi)

    ## create arguments for updateBetaH
    uniqY <- unique(y)  # y is already sorted
    yCt <- as.vector(table(y))  # yCt is the number failures at each unique failure time

    if (hazRestrict=="none") {
        coxX <- x
        coxY <- survival::Surv(y, d)
        coxInt <- NULL  # no intercept columns
    } else {
        coxX <- matrix(0, nrow=k*n, ncol=k*p)
        for (m in 1:k) coxX[((m-1)*n+1):(m*n), ((m-1)*p+1):(m*p)] <- x
        coxY <- survival::Surv(rep(y, k), rep(d, k))
        if (hazRestrict=="proportional") {
            ## intercept columns:
            coxInt <- sapply(1:k, function(m) c(rep(0, n*(m-1)), rep(1, n), rep(0, n*(k-m))))
                # class 1 is baseline for updateBetaH function
        } else {
            coxInt <- NULL  # no intercept columns
        }
    }

    ## Compute EM operator M(.)
    mPi <- colMeans(cacFitProfiled$uHat[, -k, drop=FALSE])
    mbetah <- updateBetaH(hazRestrict, coxX, coxY, coxInt, uniqY, yCt, cacFitProfiled$uHat,
                            fixedBeta=NULL, offsetMat=NULL)
    mBeta <- c(mbetah$betaHat)
    if (hazRestrict == "proportional") {
        hHat <- mbetah$hHat
        mInt <- colMeans(log(hHat / hHat[, k]), na.rm=TRUE)  # class k is baseline for inference
        mBeta <- c(mBeta, mInt[-k])
    }

    ## Compute score and observed info of profiled expected complete data loglik
    si <- cacProfileScoreInfoECD(cacFitProfiled)
    profileScorePi <- si$profileScorePi
    profileScoreBeta <- si$profileScoreBeta
    profileInfoECDPi <- si$profileInfoECDPi
    profileInfoECDBeta <- si$profileInfoECDBeta

    list(mPi=mPi, mBeta=mBeta, profileScorePi=profileScorePi, profileScoreBeta=profileScoreBeta,
         profileInfoECDPi=profileInfoECDPi, profileInfoECDBeta=profileInfoECDBeta)
}
