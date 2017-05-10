#' Random starting weights
#'
#' Constructor function for the list of parameters to be passed to the \code{uRandomControl}
#' argument of the \code{cac} function. For this starting method, a class is randomly
#' selected for each observation and given starting weight \code{p}. The remaining
#' weight is divided among the other classes.
#'
#' @param nStart Number of random starts.
#' @param p The amount of weight to be assigned to a randomly selected class
#' for each observation. The remaining weight is split between the remainning classes.
#' \code{p} should be between zero and one, and is intended to be closer to one.
#'
#' @return An S3 object of class 'uRandomControl' containing \code{nStart} and \code{p}.
#'
#' @export
getURandomControl <- function(nStart=10, p=.8) {
    list(nStart=nStart, p=p)
}

#' Random subset starting values
#'
#' Constructor function for the list fo parameters to be passed to the \code{uSubsetControl}
#' argument of the \code{cac} function. This method randomly selects a subset of
#' observations of size \code{npp} for each class. Beta for each class is initialized
#' to the Cox model fit within the corresponding subset. Observations are drawn
#' without replacement. The mixture probabilities are initialized to be equal, and
#' the baseline hazard functions are each initialized to be the
#' Nelson-Aalen estimate from the full-data.
#'
#' @param nStart Number of random starts.
#' @param npp Number of observations randomly selecting for each starting subset.
#' @param maxiter Maximum attempts before resorting to getUStartRandom. It is possible
#' for the "subset" starting method to fail, for example if one of the subsets
#' contains only censored observations. If "subset" starts faill too frequently,
#' it may help to increase \code{npp}.
#'
#' @return An S3 object of class 'uSubsetControl' containing \code{nStart},
#' \code{npp}, and \code{maxiter}.
#'
#' @export
getUSubsetControl <- function(nStart=10, npp=NULL, maxiter=100) {
    list(nStart=nStart, npp=npp, maxiter=maxiter)
}

getUStartRandom <- function(n, k, nStart, p)
{
    if (k==1) {
        uStart <- replicate(nStart, matrix(1, nrow=n, ncol=1), simplify=FALSE)
    } else {
        uStart <- replicate(nStart, simplify=FALSE, {
            startClass <- sample(1:k, n, replace=TRUE)
            uStart <- matrix( (1-p) / (k-1), nrow=n, ncol=k)
            uStart[cbind(1:n, startClass)] <- p
            uStart
        })
    }
    uStart
}

# try setting lambda and paritionSize stochastically
# hone in on values that get best objective function results
# give uncensored observations higher probability of being sampled?
getUStartSubset <- function(x, y, d, k, uThresh, npp, nStart, maxiter)
{
    n <- nrow(x)
    p <- ncol(x)
    coxY <- survival::Surv(y, d)
    piHat <- rep(1/k, k)

    # replace p with effective df for npp default calculation?
    if (is.null(npp)) npp <- min(round(3*p/mean(d)), floor(n/k))

    ## Use overall Nelson-Aalen estimator for baseline hazard in all classes,
    ## regardless of whether or not the common hazard assumption is used.
    ## Need probability mass at each failure time in each class
    ## survfit(coxph(Surv(y, d)~1)) gives Nelson-Aalen
    ## survfit(Surv(y, d)~1) gives Kaplan-Meier
    baseMod <- survival::coxph(coxY~x)
    sf <- survival::survfit(baseMod)
    timeIndex <- match(coxY[,"time"], sf$time)
    hHat <- with(sf, (n.event/n.risk)[timeIndex])
    HHat <- sf$cumhaz[timeIndex]

    uStart <- replicate(nStart, simplify=FALSE, {
        continue <- TRUE
        iter <- 0
        while (continue && iter<maxiter)
        {
            iter <- iter + 1
            continue <- tryCatch({
                ## only sample from censored observations if necessary
                # d1 <- which(d==1)
                # d0 <- sample(which(d==0), max(0, npp-sum(d==1)))
                # samp <- sample(c(d1, d0), npp*k)

                ## sample from all observations
                samp <- sample(n, k*npp)

                partID <- rep(1:k, each=npp)
                partLS <- lapply(1:k, function(i) samp[partID==i])
                betaLS <- lapply(partLS, function(ss) {
                    suppressWarnings(coef(survival::coxph(coxY~x, subset=ss)))})
                betaHat <- do.call(cbind, betaLS)
                xb <- x %*% betaHat # n x k matrix
                expXB <- exp(xb)
                uHat <- updateU(xb, expXB, hHat, HHat, d, piHat, uThresh)
                continue <- any(is.na(uHat))
                continue
            }, error = function(e) TRUE)

            if (continue)
                message(paste("uStart failure", iter))
        }

        if (continue)  # iteration limit exceeded
            uHat <- getUStartRandom(n, k, nStart=1, p=.8)[[1]]

        uHat
    })
    uStart
}
