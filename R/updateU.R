## xb only used by log-sum-exp version of updateU

# updateU <- function(xb, expXB, hHat, HHat, d, piHat, uThresh)
# {
#     n <- nrow(expXB)
#     hexpXB <- hHat * expXB
#     HexpXB <- HHat * expXB
#     uNum <- rep(piHat, each=n) * (hexpXB)^d * exp(-HexpXB)
#     uHat <- uNum / rowSums(uNum)
#     belowThreshID <- uHat < uThresh
#     if (any(belowThreshID)) {
#         uHat[belowThreshID] <- 0
#         uHat <- uHat / rowSums(uHat)
#     }
#     uHat
# }

updateU <- function(xb, expXB, hHat, HHat, d, piHat, uThresh)
{
    n <- nrow(xb)
    HexpXB <- HHat * expXB
    logTerm <- log(hHat) + xb
    logTerm[d==0] <- 0
    loguNum <- logTerm + log(rep(piHat, each=n)) - HexpXB
    loguDen <- apply(loguNum, 1, logSumExp, xLogScl=TRUE)
    uHat <- exp(loguNum - loguDen)
    belowThreshID <- uHat < uThresh
    if (any(belowThreshID)) {
        uHat[belowThreshID] <- 0
        uHat <- uHat / rowSums(uHat)
    }
    uHat
}

# if xLogScl=FALSE, x are exponential terms to sum: logSumExp(x) = log(sum(x))
# if xLogScl=TRUE, x are terms to exponentiate and sum: logSumExp(x) = log(sum(exp(x)))
logSumExp <- function(x, xLogScl=FALSE)
{
    if (!xLogScl) x <- log(x)
    i <- which.max(x)
    m <- x[i]
    if (m == -Inf) return(-Inf)
    if (m == Inf) return(Inf)
    lse <- log1p(sum(exp(x[-i]-m))) + m
    return(lse)
}
