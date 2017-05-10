# The complete log-likelihood can be written:
# log L(pi, beta, h; y, delta, u | x) = log L_1(pi; u) + log L_2(beta, h; y, delta | u, x)
# pi and (beta, h) are asymptotically independent
# Here, we profile L_2 over h and compute the information of the profiled complete data log-likelihood
# Future work: consider the full information of the complete data log-likelihood with (pi, beta, h)

louisInfoPi <- function(cacFit)
{
    k <- ncol(cacFit$uHat)
    u <- cacFit$uHat[, -k, drop=FALSE]
    uk <- cacFit$uHat[, k]
    pi <- cacFit$piHat[-k]
    pik <- cacFit$piHat[k]

    m1 <- sapply(1:(k-1), function(i) {
        v <- sapply(i:(k-1), function(j) sum(tcrossprod(u[, i], u[, j])))
        c(rep(NA, i-1), v)
    })
    m1[upper.tri(m1)] <- t(m1)[upper.tri(m1)]
    m1 <- (m1 - crossprod(u)) / tcrossprod(pi)

    m2 <- sapply(1:(k-1), function(i) sum(tcrossprod(u[, i], uk)))
    m2 <- (m2 - crossprod(u, uk)) / (pi*pik)
    m2 <- matrix(m2, nrow=k-1, ncol=k-1, byrow=TRUE)

    m3 <- (sum(tcrossprod(uk)) - crossprod(uk)) / pik^2
    m3 <- matrix(m3, nrow=k-1, ncol=k-1)

    infoPi <- m2 + t(m2) - m1 - m3
    infoPi
}

louisInfoBeta <- function(cacFit)
{
    x <- cacFit$x
    y <- cacFit$y
    d <- cacFit$d
    betaHat <- cacFit$betaHat
    uHat <- cacFit$uHat
    n <- nrow(x)
    k <- ncol(betaHat)

    siList <- lapply(1:k, function(i)
        coxScoreInfo(betaHat[, i], x, y, d, wt=uHat[, i], offset=NULL))
    infoTerm <- Reduce("+", lapply(siList, function(x) x$info))

    scoreMatList <- lapply(siList, function(x) do.call(rbind, x$scoreTerms))
    scoreOP1 <- Reduce("+", lapply(1:k, function(i)
        crossprod(uHat[, i] * scoreMatList[[i]], scoreMatList[[i]])))
    id1 <- rep(1:(n-1), times=(n-1):1)
    id2 <- do.call(c, lapply(2:n, function(i) i:n))
    scoreMat <- Reduce("+", scoreMatList)
    scoreOP2 <- crossprod(scoreMat[id1, ], scoreMat[id2, ])
    info <- infoTerm - (scoreOP1 + scoreOP2 + t(scoreOP2))

    scoreOPTerm <- Reduce("+", lapply(siList, function(x) tcrossprod(x$score)))
    info2 <- infoTerm - scoreOPTerm

    list(info=info, info2=info2)
}

louisInfo <- function(cacFit)
{
    piInfo <- louisInfoPi(cacFit)
    betaInfo <- louisInfoBeta(cacFit)$info
    piVar <- solve(piInfo)
    betaVar <- solve(betaInfo)
    piSE <- sqrt(diag(piVar))
    betaSE <- sqrt(diag(betaVar))
    list(piSE=piSE, betaSE=betaSE, piVar=piVar, betaVar=betaVar,
         piInfo=piInfo, betaInfo=betaInfo)
}
