# Function to return score and observed information of the expected complete data (ECD) log-likelihood, profiled over h

# Note: for hazRestrict=="proportional", the intercepts are the first (k-1) terms of the score and information.

# If the expectation (of uHat) is taken at the current parameter values (pi, beta, h),
# then the score of the ECD equals the score of the observed data (McLachlan eq. 3.15 and Jamshidian Jennrich eq. 2).
# At EM convergence, the score is zero.
# The ECD observed information has independent pi and beta components. I do not believe the CAC
# estimators piHat and betaHat are necessarily independent, however.
# The score of the ECD (and observed data), at points surrounding the EM solution, is used for the NDS method.
# The information of the ECD, at the EM solution, along with EM updates at surrounding points, is used for the NDM (i.e. SEM) method.

cacProfileScoreInfoECDPi <- function(cacFit)
{
    n <- nrow(cacFit$uHat)
    k <- ncol(cacFit$uHat)
    u <- cacFit$uHat[, -k, drop=FALSE]
    uk <- cacFit$uHat[, k]
    pi <- cacFit$piHat[-k]
    pik <- cacFit$piHat[k]

    # profileScore is the score of the observed data log-likelihood (profiled over h), calculated as
    # the score of the expected complete data (ECD) log-likelihood (profiled over h)
    profileScoreTerms <- u / rep(pi, each=n) - uk / pik  # n x (k-1) matrix
    profileScore <- colSums(profileScoreTerms)

    # profileInfoECD is the observed information of the ECD log-likelihood (profiled over h)
    profileInfoECD <- diag(colSums(u) / pi^2, nrow=k-1) + sum(uk) / pik^2  # this constant shows up in all terms

    list(profileScorePi=profileScore, profileInfoECDPi=profileInfoECD)
}

cacProfileScoreInfoECDBeta <- function(cacFit)
{
    n <- nrow(cacFit$uHat)
    k <- ncol(cacFit$uHat)
    p <- length(cacFit$betaHat)
    uHat <- cacFit$uHat
    betaHat <- cacFit$betaHat
    arglist <- cacFit$arglist
    x <- arglist$x
    y <- arglist$y
    d <- arglist$d
    hazRestrict <- arglist$hazRestrict

    if (hazRestrict == "none")
    {
        siList <- lapply(1:k, function(i) coxProfileScoreInfo(betaHat[, i], x, y, d, wt=uHat[, i]))
        profileScore <- do.call(c, lapply(siList, function(si) si$profileScore))
        profileInfoECD <- matrix(0, nrow=k*p, ncol=k*p)
        for (i in 1:k) profileInfoECD[((i-1)*p+1):(i*p), ((i-1)*p+1):(i*p)] <- siList[[i]]$profileInfo
    } else
    {
        coxX <- matrix(0, nrow=k*n, ncol=k*p)
        for (i in 1:k) coxX[((i-1)*n+1):(i*n), ((i-1)*p+1):(i*p)] <- x
        coxY <- rep(y, k)
        coxD <- rep(d, k)
        coxWt <- as.vector(uHat)
        beta <- as.vector(betaHat)
        if (hazRestrict == "identical")
        {
            si <- coxProfileScoreInfo(beta, coxX, coxY, coxD, coxWt)
        } else  # (hazRestrict == "proportional")
        {
            ## Method 1: Add intercepts
            coxInt <- do.call(cbind, lapply(1:(k-1), function(m) c(rep(0, n*(m-1)), rep(1, n), rep(0, n*(k-m)))))
                # Note: class k is baseline
            coxX <- cbind(coxX, coxInt)
            hHat <- cacFit$hHat
            int <- colMeans(log(hHat / hHat[, k]), na.rm=TRUE)[-k]
            beta <- c(beta, int)
            si <- coxProfileScoreInfo(beta, coxX, coxY, coxD, coxWt)

            ## Method 2: Calculate as if intercept estimates are true values
            ## Score and information, calculated as if the offsets are known. Not sure how well this works.
            # hHat <- cacFit$hHat
            # offsetVals <- colMeans(log(hHat / hHat[, 1]), na.rm=TRUE)
            # NA values are from 0 / 0
            # All non-NA ratios should be the same, but take colMeans anyway
            # offset <- rep(offsetVals, each=n)
            # si <- coxProfileScoreInfo(beta, coxX, coxY, coxD, coxWt, offset)

            ## Method 3: invert twice
            ## Another possible way to get the information? Not sure about score.
            # beta1 <- c(offsetVals[-1], beta)
            # coxX1 <- cbind(do.call(cbind, lapply(1:(k-1), function(i) c(rep(0, n*i), rep(1, n), rep(0, n*(k-i-1))))), coxX)
            # si1 <- coxProfileScoreInfo(beta1, coxX1, coxY, coxD, coxWt)
            # profileInfoECD <- solve(solve(si1$profileInfo)[k:(k-1+k*p), k:(k-1+k*p)])
        }

        profileScore <- si$profileScore
        profileInfoECD <- si$profileInfo
    }

    list(profileScoreBeta=profileScore, profileInfoECDBeta=profileInfoECD)
}

cacProfileScoreInfoECD <- function(cacFit)
{
    siPi <- cacProfileScoreInfoECDPi(cacFit)
    siBeta <- cacProfileScoreInfoECDBeta(cacFit)
    c(siPi, siBeta)
}
