## Score and observed information using Breslow likelihood for ties
## offset can be either a single number or length n vector
coxProfileScoreInfo <- function(beta, x, y, d, wt=NULL, offset=NULL)
{
    if (!is.matrix(x)) stop("x must be a matrix")

    n <- nrow(x)
    p <- ncol(x)
    if (is.null(wt)) wt <- rep(1, n)
    if (is.null(offset)) offset <- rep(0, n)

    if (length(y) != n) stop("length(y) != nrow(x)")
    if (length(d) != n) stop("length(d) != nrow(x)")
    if (length(beta) != p) stop("length(beta) != ncol(x)")
    if (length(wt) != n) stop("length(wt) != nrow(x)")
    if (length(offset) != n) stop("length(offset) != nrow(x)")

    hr <- exp(c(x %*% beta) + offset)

    siTerms <- lapply(1:n, function(i) {

        if (d[i]==0 || wt[i]==0) {
            scoreTerm <- rep(0, p)
            infoTerm <- matrix(0, nrow=p, ncol=p)
            return(list(scoreTerm=scoreTerm, infoTerm=infoTerm))
        }

        riskSet <- which(y>=y[i])
        nRisk <- length(riskSet)
        xRisk <- x[riskSet, , drop=FALSE]
        hrRisk <- wt[riskSet] * hr[riskSet]
        hazWt <- hrRisk / sum(hrRisk)
            # weight of each observation at risk at failure time t
            # weights are conditional failure probability for each observation

        xRiskWtd <- hazWt * x[riskSet, , drop=FALSE]
        wtdAvgX <- colSums(xRiskWtd)
            # wtd avg of x vectors in risk set
        wtdAvgXX <- crossprod(xRisk, xRiskWtd)
            # wtd avg of xx' matrices in risk set

        scoreTerm <- wt[i] * (x[i, ] - wtdAvgX)
        infoTerm <- wt[i] * (wtdAvgXX - tcrossprod(wtdAvgX))

        list(scoreTerm=scoreTerm, infoTerm=infoTerm)
    })

    profileScoreTerms <- lapply(siTerms, function(x) x$scoreTerm)
    profileInfoTerms <- lapply(siTerms, function(x) x$infoTerm)
    profileScore <- Reduce("+", profileScoreTerms)
    profileInfo <- Reduce("+", profileInfoTerms)
    list(profileScore=profileScore, profileInfo=profileInfo,
         profileScoreTerms=profileScoreTerms, profileInfoTerms=profileInfoTerms)
}
