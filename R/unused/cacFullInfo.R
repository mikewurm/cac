# Not completed. This was an attempt to calculate the full (non-profile) information for CAC.
# There seems to be an error or some numerical instability because the information is not always positive definite.
# The information also usually has a high condition number.

# infoECD is the observed info (negative Hessian) of the expected complete data, which can be calculated by summing over all observations
# scoreOP objects are the sum over all observations of the expected outer product of the complete data score
# tcrossprod(scoreTerms) are the sum over all observations of the expected complete data score

cacFullInfo <- function(cacFit)
{
    arglist <- cacFit$arglist
    hazRestrict <- arglist$hazRestrict
    x <- arglist$x
    y <- arglist$y
    d <- arglist$d

    piHat <- cacFit$piHat
    betaHat <- cacFit$betaHat
    hHat <- cacFit$hHat
    HHat <- cacFit$HHat
    hHatOrd <- cacFit$hHatOrd
    HHatOrd <- cacFit$HHatOrd
    uHat <- cacFit$uHat

    n <- nrow(x)
    p <- ncol(x)
    k <- arglist$k
    u <- cacFit$uHat[, -k, drop=FALSE]
    uk <- cacFit$uHat[, k]
    pi <- cacFit$piHat[-k]
    pik <- cacFit$piHat[k]
    hr <- exp(x %*% betaHat)  # n x k matrix
    yy <- sort(unique(y[d==1]))

    # Pi block
    scoreTermsPi <- u / rep(pi, each=n) - uk / pik  # n x (k-1) matrix
    infoPi <- crossprod(scoreTermsPi)
    # Alternate info calculation:
    # infoECDPi <- Reduce(`+`, lapply(1:n, function(i) diag(u[i, ] / pi^2) + uk[i] / pik^2))
    # scoreOPPi <- Reduce(`+`, lapply(1:n, function(i) diag(u[i, ] / pi^2) + uk[i] / pik^2))
    # infoPi <- infoECDPi - scoreOPPi + crossprod(scoreTermsPi)

    # Beta block
    # scoreOPBeta is a k-block diagonal matrix
    # t1 <- uHat * HHat * hr
    # t1 <- HHat * hr
    # t2 <- d - HHat * hr
    scoreTermsBeta <- do.call(cbind, lapply(1:k, function(m) x * (uHat[, m] * (d - HHat[, m] * hr[, m]))))
    infoECDBeta <- matrix(0, nrow=k*p, ncol=k*p)
    scoreOPBeta <- matrix(0, nrow=k*p, ncol=k*p)
    for (m in 1:k) {
        index <- 1:p + p * (m-1)
        # infoECDBeta[index, index] <- crossprod(uHat[, m] * t1[, m] * x, t1[, m] * x)
        # scoreOPBeta[index, index] <- crossprod(uHat[, m] * t2[, m] * x, t2[, m] * x)
        infoECDBeta[index, index] <- Reduce(`+`, lapply(1:n, function(i) {
            uHat[i, m] * HHat[i, m] * hr[i, m] * tcrossprod(x[i, ])
        }))
        scoreOPBeta[index, index] <- Reduce(`+`, lapply(1:n, function(i) {
            uHat[i, m] * (d[i] - HHat[i, m] * hr[i, m])^2 * tcrossprod(x[i, ])
        }))
    }
    infoBeta <- infoECDBeta - scoreOPBeta + crossprod(scoreTermsBeta)

    # Pi-beta block:
    # infoECD is a matrix of zeros
    scoreOPPiBeta <- matrix(0, nrow=k-1, ncol=k*p)
    indexk <- 1:p + p*(k-1)
    for (m in 1:(k-1)) {
        index <- 1:p + p * (m-1)
        # scoreOPPiBeta[i, index] <- c(crossprod(x, uHat[, m] * t2[, m] / pi[m]))
        # scoreOPPiBeta[m, indexk] <- -c(crossprod(x, uHat[, k] * t2[, k] / pik))
        scoreOPPiBeta[m, index] <- Reduce(`+`, lapply(1:n, function(i) {
            x[i, ] * uHat[i, m] * (d[i] - HHat[i, m] * hr[i, m]) / pi[m]
        }))
        scoreOPPiBeta[m, indexk] <- Reduce(`+`, lapply(1:n, function(i) {
            x[i, ] * uHat[i, k] * (d[i] - HHat[i, k] * hr[i, k]) / pik
        }))
    }
    infoPiBeta <- -scoreOPPiBeta + crossprod(scoreTermsPi, scoreTermsBeta)

    if (hazRestrict == "none")
    {
        # h block:
        hStar <- hHat[match(yy, y), ]  # hazard estimate at each unique failure time
        scoreTermsh <- do.call(cbind, lapply(1:k, function(m) {
            do.call(rbind, lapply(1:n, function(i) {
                uHat[i, m] * (y[i]==yy) * (d[i]==1) / (hStar[, m]+1e-100) - (y[i]>=yy) * hr[i, m] * uHat[i, m]
            }))
        }))
        infoECDh <- diag(do.call(c, lapply(1:k, function(m) {
            Reduce(`+`, lapply(1:n, function(i) {
                uHat[i, m] * (y[i]==yy) * (d[i]==1) / (hStar[, m]^2+1e-100)
            }))
        })))
        scoreOPh <- matrix(0, nrow=length(hStar), ncol=length(hStar))
        for (m in 1:k) {
            index <- 1:nrow(hStar) + nrow(hStar) * (m-1)
            scoreOPh[index, index] <- Reduce(`+`, lapply(1:n, function(i) {
                diag(uHat[i, m] * (y[i] == yy) * (d[i] == 1) / (hStar[, m]^2+1e-100)) -  # cancels with infoECDh
                    2 * diag(uHat[i, m] * (y[i] == yy) * (d[i] == 1) / (hStar[, m]+1e-100) * hr[i, m]) +
                    tcrossprod(y[i] >= yy) * uHat[i, m] * hr[i, m]^2
            }))
        }
        infoh <- infoECDh - scoreOPh + crossprod(scoreTermsh)

        # Pi-h block (k-1) x k*D
        # infoECDPih is a matrix of zeros
        scoreOPPih <- matrix(0, nrow=k-1, ncol=length(hStar))
        indexk <- 1:nrow(hStar) + nrow(hStar) * (k-1)
        for (m in 1:(k-1)) {
            index <- 1:nrow(hStar) + nrow(hStar) * (m-1)
            scoreOPPih[m, index] <- Reduce(`+`, lapply(1:n, function(i) {
                (y[i] == yy) * (d[i] == 1) / (hStar[, m]+1e-100) * (u[i, m] / pi[m]) -
                    (y[i] >= yy) * (u[i, m] / pi[m] * hr[i, m])
            }))
            scoreOPPih[m, indexk] <- Reduce(`+`, lapply(1:n, function(i) {
                (y[i] == yy) * (d[i] == 1) / (hStar[, m]+1e-100) * (-uk[i] / pik) -
                    (y[i] >= yy) * (-uk[i] / pik * hr[i, k])
            }))
        }
        infoPih <- -scoreOPPih + crossprod(scoreTermsPi, scoreTermsh)

        # Beta-h block k*p x k*D
        infoECDBetah <- scoreOPBetah <- matrix(0, nrow=k*p, ncol=length(hStar))
        for (m in 1:k) {
            index1 <- 1:p + p * (m-1)
            index2 <- 1:nrow(hStar) + nrow(hStar) * (m-1)
            infoECDBetah[index1, index2] <- Reduce(`+`, lapply(1:n, function(i) {
                tcrossprod(x[i, ], (y[i] >= yy)) * uHat[i, m] * hr[i, m]
            }))
            scoreOPBetah[index1, index2] <- Reduce(`+`, lapply(1:n, function(i) {
                tcrossprod(x[i, ] * uHat[i, m] * (d[i] - HHat[i, m] * hr[i, m]),
                           (y[i] == yy) * (d[i] == 1) / (hStar[, m]+1e-100) - (y[i] >= yy) * hr[i, m])
            }))
        }
        infoBetah <- infoECDBetah - scoreOPBetah + crossprod(scoreTermsBeta, scoreTermsh)
    } else  # (hazRestrict %in% c("identical", "proportional")
    {
        # h block:
        nFail <- sapply(yy, function(t) sum(y == t))  # number of failures at each unique failure time
        hStar <- hHat[match(yy, y), 1]  # hazard estimate at each unique failure time
        scoreTermsh <- do.call(rbind, lapply(1:n, function(i) (y[i]==yy) * (d[i]==1) / hStar - (y[i]>=yy) * sum(hr[i, ] * uHat[i, ])))
        infoECDh <- diag(nFail / hStar^2)
        scoreOPh <- Reduce(`+`, lapply(1:n, function(i) {
            diag((y[i] == yy) * (d[i] == 1) / hStar^2) -  # cancels with infoECDh
                2 * diag((y[i] == yy) * (d[i] == 1) / hStar * sum(uHat[i, ] * hr[i, ])) +
                tcrossprod(y[i] >= yy) * sum(uHat[i, ] * hr[i, ]^2)
        }))
        infoh <- infoECDh - scoreOPh + crossprod(scoreTermsh)

        # Pi-h block (k-1) x D
        # infoECDPih is a matrix of zeros
        scoreOPPih <- do.call(rbind, lapply(1:(k-1), function(m) {
            Reduce(`+`, lapply(1:n, function(i) {
                (y[i] == yy) * (d[i] == 1) / hStar * (u[i, m] / pi[m] - uk[i] / pik) -
                    (y[i] >= yy) * (u[i, m] / pi[m] * hr[i, m] - uk[i] / pik * hr[i, k])
            }))
        }))
        infoPih <- -scoreOPPih + crossprod(scoreTermsPi, scoreTermsh)

        # Beta-h block k*p x D
        infoECDBetah <- do.call(rbind, lapply(1:k, function(m) {
            Reduce(`+`, lapply(1:n, function(i) {
                tcrossprod(x[i, ], (y[i] >= yy)) * uHat[i, m] * hr[i, m]
            }))
        }))
        scoreOPBetah <- do.call(rbind, lapply(1:k, function(m) {
            Reduce(`+`, lapply(1:n, function(i) {
                tcrossprod(x[i, ] * uHat[i, m] * (d[i] - HHat[i, m] * hr[i, m]),
                           (y[i] == yy) * (d[i] == 1) / hStar - (y[i] >= yy) * hr[i, m])
            }))
        }))
        infoBetah <- infoECDBetah - scoreOPBetah + crossprod(scoreTermsBeta, scoreTermsh)

        if (hazRestrict == "proportional")
        {
            # zeta block

            # zeta-pi block

            # zeta-beta block

            # zeta-h block
        }
    }

    info <- rbind(cbind(infoPi, infoPiBeta, infoPih),
                  cbind(t(infoPiBeta), infoBeta, infoBetah),
                  cbind(t(infoPih), t(infoBetah), infoh))
    hDim <- ncol(info) - (k-1) - k*p
    lincom1 <- diag(k-1+k*p+hDim)
    lincom <- rbind(lincom1[1:(k-1), ], rep(c(-1, 0), c(k-1, k*p+hDim)), lincom1[k:(k-1+k*p+hDim), ])
    varAll <- lincom %*% solve(info) %*% t(lincom)
    varPi <- varAll[1:k, 1:k]
    varBeta <- varAll[(k+1):(k+k*p), (k+1):(k+k*p)]
    varh <- varAll[(k+k*p+1):(k+k*p+hDim), (k+k*p+1):(k+k*p+hDim)]
    sePi <- sqrt(diag(varPi))
    seBeta <- matrix(sqrt(diag(varBeta)), ncol=k)

    list(sePi=sePi, seBeta=seBeta,
         varPi=varPi, varBeta=varBeta, varh=varh, varAll=varAll, info=info)
}

