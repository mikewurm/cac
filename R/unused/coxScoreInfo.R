# Does not calculate individual score and info terms needed for Louis formula

## Score and observed information using Breslow likelihood for ties
coxScoreInfo <- function(beta, x, y, d, wt=NULL, offset=NULL)
{
    if (is.null(wt)) wt <- rep(1, nrow(x))
    if (is.null(offset)) offset <- rep(0, nrow(x))
    hr <- exp(as.vector(x %*% beta) + offset)
    failTimes <- unique(y[d==1 & wt>0])
    failWt <- sapply(failTimes, function(t) sum(wt[y==t]))
    nFail <- length(failTimes)

    wtdAvgTerms <- lapply(failTimes, function(t) {

        riskID <- which(y>=t)
        nRisk <- length(riskID)
        xRisk <- x[riskID, , drop=FALSE]
        hrRisk <- wt[riskID] * hr[riskID]
        hazWt <- hrRisk / sum(hrRisk)
        # weight of each observation at risk at failure time t
        # weights are conditional failure probability for each observation
        xRiskWtd <- hazWt * x[riskID, , drop=FALSE]

        wtdAvgX <- colSums(xRiskWtd)
        # wtd avg of x vectors in risk set
        wtdAvgXX <- Reduce("+", lapply(1:nRisk, function(i) {
            tcrossprod(xRisk[i, ], xRiskWtd[i, ])}))
        # wtd avg of xx' matrices in risk set
        list(wtdAvgX=wtdAvgX, wtdAvgXX=wtdAvgXX)
    })

    wtdAvgWtdAvgX <- Reduce("+", lapply(1:nFail, function(i) {
        failWt[i] * wtdAvgTerms[[i]]$wtdAvgX}))
    wtdAvgWtdAvgXX <- Reduce("+", lapply(1:nFail, function(i) {
        failWt[i] * wtdAvgTerms[[i]]$wtdAvgXX}))
    wtdAvgOPWtdAvgX <- Reduce("+", lapply(1:nFail, function(i) {
        failWt[i] * tcrossprod(wtdAvgTerms[[i]]$wtdAvgX)}))
    xFailSum <- colSums(wt[d==1] * x[d==1, , drop=FALSE])
    score <- xFailSum - wtdAvgWtdAvgX
    info <- wtdAvgWtdAvgXX - wtdAvgOPWtdAvgX
    list(score=score, info=info)
}
