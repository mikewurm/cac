## Simulate single PH mixture data set
#' @export
simDataPHM <- function(nObs, beta, xDist, survDist, censDist, classDist)
{
    p <- nrow(beta)
    k <- ncol(beta)

    ## Simulate data
    x <- xDist(nObs, p)
    u <- classDist(nObs, k)
    xb <- x %*% beta
    xBeta <- xb[cbind(1:nObs, u)]
    hr <- exp(xBeta)
    t <- survDist(nObs, hr)
    cc <- censDist(nObs)
    y <- pmin(cc, t)
    d <- as.integer(t < cc)
    list(u=u, x=x, y=y, d=d)
}

## Get summary metric for all scenarios
## (cacSim returns a list of metrics by scenario)
#' @export
aggSim <- function(sim, metric, digits=NULL)
{
    if (is.null(digits))
    {
        lapply(sim, function(s) s[[metric]])
    } else
    {
        lapply(sim, function(s) round(s[[metric]], digits=digits))
    }
}

## TESTING PERCENTAGE CENSORED
## could improve by getting mean/sd over multiple reps
## could also get censoring pct by class (should be the same for symmetric x & beta though)
## startVals parameter does not affect censoring pct
## nRep doesn't have to be the same as number of simulations
## Returns a matrix of percentages: rows correspond to simulation settings, columns to classes
#' @export
testPctCens <- function(nRep, nObs, beta, xDist, survDist, censDist, classDist)
{
    paramGrid <- expand.grid(nObs = seq_along(nObs),
                             beta = seq_along(beta),
                             xDist = seq_along(xDist),
                             survDist = seq_along(survDist),
                             censDist = seq_along(censDist),
                             classDist = seq_along(classDist))

    pctCens <- apply(paramGrid, 1, function(par)
    {
        pctCensReps <- replicate(nRep, {
            nObs <- nObs[[par["nObs"]]]
            beta<- beta[[par["beta"]]]
            xDist <- xDist[[par["xDist"]]]
            survDist <- survDist[[par["survDist"]]]
            censDist <- censDist[[par["censDist"]]]
            classDist <- classDist[[par["classDist"]]]

            p <- nrow(beta)
            k <- ncol(beta)

            ## Simulate data
            x <- xDist(nObs, p)
            u <- classDist(nObs, k)
            xb <- x %*% beta
            xBeta <- xb[cbind(1:nObs, u)]
            hr <- exp(xBeta)
            t <- survDist(nObs, hr)
            cc <- censDist(nObs)
            y <- pmin(cc, t)
            d <- as.integer(t < cc)
            pctCens <- sapply(1:k, function(k) 1 - mean(d[u==k]))
            pctCens
        })

        pctCens <- rowMeans(pctCensReps)
        pctCens
    })

    ## pctCens could be a list if #classes is different by scenario
    if (is.matrix(pctCens)) pctCens <- t(pctCens)
    pctCens
}
