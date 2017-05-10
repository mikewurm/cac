# This function calculates the non-profile score and information for a Cox PH model.
# It can be used for inference on both beta and baseline hazard, and the beta
# portion of the inverse information is identical to that of the profile likelihood.
# According to Klein/Moeschberger, cumulative baseline hazard standard errors are
# identical between this method and the typical Nelson-Aalen approach.
#' @export
coxFullScoreInfo <- function(beta, x, y, d)
{
    n <- length(y)
    yy <- sort(unique(y[d==1]))  # unique failure times
    nFail <- sapply(yy, function(t) sum(y == t))
    riskSetList <- lapply(yy, function(t) which(y >= t))
    hr <- exp(c(x %*% beta))
    hStar <- nFail / sapply(riskSetList, function(r) sum(hr[r]))
    HStar <- cumsum(hStar)
    h <- sapply(y, function(t) sum(hStar[yy == t]))
    H <- sapply(y, function(t) sum(hStar[yy <= t]))

    # Score
    # scoreTermsBeta <- x * (d - H * hr)
    # scoreBeta <- colSums(scoreTermsBeta)
    scoreBeta <- c(crossprod(x, d - H * hr)) # alternate
    scoreTermsh <- t(sapply(1:n, function(i) (y[i]==yy && d[i]==1) / hStar - hr[i] * (y[i] >= yy)))
    scoreh <- colSums(scoreTermsh)
    fullScore <- c(scoreBeta, scoreh)

    # Info
    # infoTermsBeta <- lapply(1:n, function(i) H[i] * hr[i] * tcrossprod(x[i, ]))
    # infoBeta <- Reduce("+", infoTermsBeta)
    infoBeta <- crossprod(H * hr * x, x)
    # infoTermsBetah <- lapply(1:n, function(i) hr[i] * tcrossprod(x[i, ], (y[i] >= yy)))
    # infoBetah <- Reduce("+", infoTermsBetah)
    infoBetah <- crossprod(hr * x, outer(y, yy, ">="))
    # infoTermsh <- lapply(1:n, function(i) diag(d[i] / hStar^2 * (y[i] == yy)))
    # infoh <- Reduce("+", infoTermsh)
    infoh <- diag(nFail / hStar^2)
    fullInfo <- rbind(cbind(infoBeta, infoBetah), cbind(t(infoBetah), infoh))

    # Test score correctness
    # loglik <- function(beta) sum(log(h[d==1])) + sum(x[d==1, ] %*% beta) - sum(H * exp(x %*% beta))
    # loglik(beta)
    # eps <- 1e-6
    # del <- c(1, rep(0, length(beta)-1)) * eps
    # (loglik(beta + del) - loglik(beta)) / eps

    list(fullScore=fullScore, fullInfo=fullInfo)
}
