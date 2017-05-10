cacSingleStart <- function(uStart,
                           n, p, coxX, coxY, coxInt, uniqY, yCt,
                           x, y, d, k,
                           hazRestrict=c("none", "proportional", "identical"),
                           maxiter=Inf, eps=1e-10, uThresh=1e-8, trace=FALSE,
                           fixedBeta=NULL, fixedPi=NULL)
{
    if (!is.null(fixedPi))
        piHat <- fixedPi

    if (!is.null(fixedBeta)) {
        if (hazRestrict == "proportional") {
            intHat <- fixedBeta[1, ]
            betaHat <- fixedBeta[-1, , drop=FALSE]
            xb <- x %*% betaHat
            expXB <- exp(xb)
            offsetMat <- xb + matrix(intHat, nrow=n, ncol=k, byrow=TRUE)
        } else {
            betaHat <- fixedBeta
            xb <- x %*% betaHat
            expXB <- exp(xb)
            offsetMat <- xb
        }
    } else {
        offsetMat <- NULL
    }

    uHat <- uStart
    loglik <- -Inf
    conv <- FALSE
    iter <- 0
    while(!conv && iter<maxiter)
    {
        iter <- iter + 1

        ## E-step
        if (iter > 1) {
            uHat <- updateU(xb, expXB, hHat, HHat, d, piHat, uThresh)
        }

        ## M-step
        if (is.null(fixedPi)) {
            piHat <- colMeans(uHat)
        }

        betaH <- updateBetaH(hazRestrict, coxX, coxY, coxInt, uniqY, yCt, uHat,
                             fixedBeta, offsetMat)
        hHat <- betaH$hHat
        HHat <- betaH$HHat
        if (is.null(fixedBeta)) {
            betaHat <- betaH$betaHat
            xb <- x %*% betaHat
            expXB <- exp(xb)
        }

        if (min(expXB)==0 || max(expXB)==Inf) {
            loglik <- -Inf
            break
        }

        ## Check convergence
        # Observed data log likelihood:
        loglikOld <- loglik
        loglik <- sum(log(rowSums(rep(piHat, each=n)*exp(-expXB*HHat) * (expXB*hHat)^d)))
        dif <- abs((loglik - loglikOld) / (loglik + 1e-100))
        conv <- dif < eps

        if (trace) cat("iter ", iter, ": loglik relative change = ", dif, "\n", sep='')

    }

    # if (!conv) return(NULL)

    # Expected complete data log likelihood:
    # ll1 <- sum(uHat * rep(log(piHat), each=n))
        # careful for log(hHat)=-Inf and d=0 or u=0
    # ll2 <- sum(d * uHat * log(hHat), d * uHat * xb, -HHat * uHat * expXB, na.rm=TRUE)
    # loglikEC <- ll1 + ll2

    cacLocal <- list(conv=conv, iter=iter, betaHat=betaHat, piHat=piHat,
                     uHat=uHat, hHat=hHat, HHat=HHat, loglik=loglik)
    class(cacLocal) <- "cacLocal"
    cacLocal
}
