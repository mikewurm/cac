# Alternate version: assumes (pi, beta) are independent
# When hazRestrict = "none", also assumes that beta blocks are independent
cacProfileInfo <- function(cacFit, eta=1e-4)
{
    betaHat <- cacFit$betaHat
    piHat <- cacFit$piHat
    p <- nrow(betaHat)
    k <- ncol(betaHat)
    hazRestrict <- cacFit$arglist$hazRestrict

    reFact <- c(-2, -1, 1, 2)  # Richardson extrapolation factors
    dmPi <- dscorePi <- matrix(nrow=k-1, ncol=k-1)
    dmBeta <- dscoreBeta <- matrix(nrow=k*p, ncol=k*p)
    mListPi <- scoreListPi <- rep(list(matrix(nrow=k-1, ncol=k-1)), 4)
    mListBeta <- scoreListBeta <- rep(list(matrix(nrow=k*p, ncol=k*p)), 4)
    hPi <- rep(NA, k-1)
    hBeta <- rep(NA, k*p)

    # Calculate M(pi) and Score(pi) at Richardson Extrapolation points surrounding piHat
    for (i in 1:(k-1))
    {
        delPi <- rep(0, k)
        delPi[i] <- hPi[i] <- eta
        delPi[k] <- -hPi[i]

        for (q in 1:4) {
            msiTemp <- cacProfileMSI(beta = betaHat,
                                     pi = piHat + reFact[q] * delPi,
                                     cacFit = cacFit)
            mListPi[[q]][, i] <- msiTemp$mPi
            scoreListPi[[q]][, i] <- msiTemp$profileScorePi
        }
    }

    # Calculate M(beta) and Score(beta) at Richardson Extrapolation points surrounding betaHat
    for (i in 1:(k*p))
    {
        delBeta <- rep(0, k*p)
        delBeta[i] <- hBeta[i] <- max(abs(betaHat[i]), 1) * eta

        for (q in 1:4) {
            msiTemp <- cacProfileMSI(beta = betaHat + reFact[q] * delBeta,
                                     pi = piHat,
                                     cacFit = cacFit)
            mListBeta[[q]][, i] <- c(msiTemp$mBeta)
            scoreListBeta[[q]][, i] <- c(msiTemp$profileScoreBeta)
        }
    }

    dmPi <- (mListPi[[1]] - 8*mListPi[[2]] + 8*mListPi[[3]] - mListPi[[4]]) / (12*rep(hPi, each=k-1))
    dscorePi <- (scoreListPi[[1]] - 8*scoreListPi[[2]] + 8*scoreListPi[[3]] - scoreListPi[[4]]) / (12*rep(hPi, each=k-1))

    dmBeta <- (mListBeta[[1]] - 8*mListBeta[[2]] + 8*mListBeta[[3]] - mListBeta[[4]]) / (12*rep(hBeta, each=k*p))
    dscoreBeta <- (scoreListBeta[[1]] - 8*scoreListBeta[[2]] + 8*scoreListBeta[[3]] - scoreListBeta[[4]]) / (12*rep(hBeta, each=k*p))

    # With the non-proportional hazards assumption, the beta groups are independent
    # Set off-block-diagonal elements to zero
    if (hazRestrict == "none") {
        dmBetaTemp <- dmBeta
        dscoreBetaTemp <- dscoreBeta
        dmBeta <- dscoreBeta <- matrix(0, nrow=k*p, ncol=k*p)
        for (i in 1:k) {
            index <- 1:p + (i-1) * p
            dmBeta[index, index] <- dmBetaTemp[index, index]
            dscoreBeta[index, index] <- dscoreBetaTemp[index, index]
        }
    }

    # Get expected complete data profile information for NDM approximation,
    # as well as exact profile info (pi only)
    si <- cacProfileScoreInfo(cacFit)

    # Exact profile information (pi only)
    profileInfoPi.Exact <- si$profileInfoPi
    piLincom <- rbind(diag(k-1), rep(-1, k-1))
    # linear combination matrix to convert (k-1 x k-1) variance matrix to (k x k)
    varPi.Exact <- piLincom %*% solve(profileInfoPi.Exact) %*% t(piLincom)
    sePi.Exact <- sqrt(diag(varPi.Exact))

    # NDM and NDS profile info approximations
    profileInfoECDPi <- si$profileInfoECDPi
    profileInfoECDBeta <- si$profileInfoECDBeta
    profileInfoPi.NDM <- profileInfoECDPi %*% (diag(k-1) - dmPi)
    profileInfoBeta.NDM <- profileInfoECDBeta %*% (diag(k*p) - dmBeta)
    profileInfoPi.NDS <- -dscorePi
    profileInfoBeta.NDS <- -dscoreBeta
    varPi.NDM <- piLincom %*% solve(profileInfoPi.NDM) %*% t(piLincom)
    varBeta.NDM <- solve(profileInfoBeta.NDM)
    varPi.NDS <- piLincom %*% solve(profileInfoPi.NDS) %*% t(piLincom)
    varBeta.NDS <- solve(profileInfoBeta.NDS)
    sePi.NDM <- sqrt(diag(varPi.NDM))
    seBeta.NDM <- sqrt(diag(varBeta.NDM))
    sePi.NDS <- sqrt(diag(varPi.NDS))
    seBeta.NDS <- sqrt(diag(varBeta.NDS))

    list(sePi.Exact=sePi.Exact, sePi.NDM=sePi.NDM, sePi.NDS=sePi.NDS,
         seBeta.NDM=seBeta.NDM, seBeta.NDS=seBeta.NDS,
         varPi.Exact=varPi.Exact, varPi.NDM=varPi.NDM, varPi.NDS=varPi.NDS,
         varBeta.NDM=varBeta.NDM, varBeta.NDS=varBeta.NDS)
}
