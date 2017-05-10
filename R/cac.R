## MAIN FUNCTION

#' Cox-assisted clustering
#'
#' Fits a Cox-assisted clustering (CAC) model.
#'
#' @importFrom survival coxph basehaz
#' @importFrom utils head tail
#' @importFrom stats coef
#'
#' @param x Covariate matrix.
#' @param y Response variable. The follow-up time of each observation.
#' @param d Censoring indicator for each observation. One to indicate failure
#' and zero to indicate censoring.
#' @param k Integer. Number of classes in the mixture model.
#' @param hazRestrict Type of restriction on the baseline hazard functions.
#' "none" allows a different nonparametric estimate of the baseline hazard in each class.
#' "proprotional" restricts the functions to be proportional, and "identical" restricts
#' them to be the same.
#' @param maxiter Maximum number of EM algorithm iterations.
#' @param eps Convergence threshold. The algorithm terminates when the relative
#' change in the semiparametric log-likelihood between successive iterations
#' drops below \code{eps}.
#' @param trace Logical. If \code{TRUE}, EM algorithm progress is printed to the terminal.
#' @param uStart An optional list of starting weight values. Each component of the list
#' should be an \eqn{n \times k} matrix, where \eqn{n} is the  number of observations
#' and \eqn{k} is the number of classes. Weights should be between zero and one, and
#' each row should sum to one.
#' @param uStartMethod Only used if \code{uStart = NULL}. This argument indicates
#' the method for generating random starting values for the EM algorithm.
#' The "subset" method initializes the coefficients by fitting a
#' Cox model to randomly selected subsets. The mixture probabilities are initialized
#' to be equal, and the baseline hazard functions are each initialized to be the
#' Nelson-Aalen estimate from the full-data. Parameters for this method are passed
#' to the \code{uSubsetControl} argument. The "random" method randomly selects
#' one class to have most of the weight for each observation. The remaining weight
#' is equally divided among the remaining observation. Parameters for this method
#' are passed to the \code{uRandomControl} argument.
#' @param uRandomControl A list of arguments used if \code{uStart = NULL} and
#' \code{uStartMethod = "random"}. Includes the number of random starts and the
#' amount of weight to be placed in the majority class for each observation.
#' This argument should be constructed with the \code{uRandomControl} function.
#' @param uSubsetControl A list of arguments used if \code{uStart = NULL} and
#' \code{uStartMethod = "subset"}. Includes the number of random starts, the number
#' of observations in each starting subset, and the maximum number of tries before
#' resorting to the "random" start method.
#' This argument should be constructed with the \code{uSubsetControl} function.
#' @param keepLocal Logical. If TRUE, the returned 'cacFit' object will include
#' a list of fitted models called 'localFits'. Each fit corresponds to a different
#' set of starting values.
#' @param fixedBeta Either NULL or a matrix of dimension p x k.
#' If not NULL, the observed semiparametric likelihood is
#' optimized holding beta fixed. Mostly intended to be used by the \code{cacProfileInfo}
#' function, but this functionality is also provided for the user.
#' @param fixedPi Either NULL or a vector of k probabilities summing to one.
#' If not NULL, the observed semiparametric likelihood is
#' optimized holding pi fixed. Mostly intended to be used by the \code{cacProfileInfo}
#' function, but this functionality is also provided for the user.
#'
#' @return S3 object of class 'cacFit', which contains information about the best
#' solution found among all sets of starting values. Also contains information about
#' the local optima found by the different sets of starting values. Note that observations
#' are sorted by failure time before optimization, so \code{uHat}, \code{hHat},
#' and \code{HHat} may correspond to their order in the \code{x}, \code{y}, and \code{d}
#' arguments. The returned \code{arglist} contains \code{x}, \code{y}, and \code{d}
#' in sorted order. The 'cacFit' object consists of the following.
#' \describe{
#'   \item{conv}{Logical indicator for whether the algorithm converged within the
#'   iteration limit for the best model fit.}
#'   \item{iter}{Number of iterations for the best model fit.}
#'   \item{betaHat}{Matrix of coefficient estimates.}
#'   \item{piHat}{Vector of mixture probability estimates.}
#'   \item{uHat}{Matrix of posterior class probabilities for each observation time.}
#'   \item{hHat}{Matrix containing the nonparametric baseline hazard estimate for
#'   each observation. This is an \eqn{n \times k} matrix, and it will contain zeros
#'   for observations that are censored.}
#'   \item{HHat}{Matrix containing the cumulative hazard estimate at each observation
#'   time. This is an \eqn{n \times k} matrix.}
#'   \item{arglist}{List of arguments passed to the \code{cac} function call.}
#'   \item{bestLocalIndex}{Index number of the best fit within the list of \code{uStart}
#'   values. Note that if starting values were generated randomly, then the list
#'   of random starting weights will be returned in \code{arglist$uStart}.}
#'   \item{loglikLocal}{Vector of semiparametric log-likelihood values for each set
#'   of starting values.}
#'   \item{localFits}{Only returned if \code{keepLocal = TRUE}. A list of parameter
#'   values for each set of starting values.}
#' }
#'
#' @export
cac <- function(x, y, d, k,
                hazRestrict=c("none", "proportional", "identical"),
                maxiter=5000, eps=1e-8, trace=FALSE,
                uStart=NULL,
                uStartMethod=c("subset", "random"),
                uRandomControl=getURandomControl(),
                uSubsetControl=getUSubsetControl(),
                keepLocal=TRUE, fixedBeta=NULL, fixedPi=NULL)
{
    ###########################################################################
    ## argument checks/mods
    ###########################################################################
    hazRestrict <- match.arg(hazRestrict)
    if (k==1 && hazRestrict=="proportional")
        stop("When k=1, use hazRestrict=c(\"none\", \"identical\"")
    uStartMethod <- match.arg(uStartMethod)
    if (!is.null(uStart) && typeof(uStart)!="list") stop("uStart should be a list.")
    x <- as.matrix(x)
    d <- as.logical(d)
    n <- nrow(x)
    p <- ncol(x)

    ## order failure times
    yOrder <- order(y)
    y <- y[yOrder]
    x <- x[yOrder, , drop=FALSE]
    d <- d[yOrder]
    if (!is.null(uStart)) uStart <- lapply(uStart, function(u) u[yOrder, , drop=FALSE])
    uniqY <- unique(y)  # y is already sorted
    yCt <- as.vector(table(y))
    ## coxph returns baseline hazard for unique times (censored & uncensored)
    ## yCt is the number of times to repeat each value
    ###########################################################################
    ## argument checks/mods
    ###########################################################################

    ###########################################################################
    ## set uStart control parameters
    ###########################################################################
    if (missing(uRandomControl)) uRandomControl <- list()
    uRandomControl <- do.call(getURandomControl, uRandomControl)
    if (missing(uSubsetControl)) uSubsetControl <- list()
    uSubsetControl <- do.call(getUSubsetControl, uSubsetControl)
    ###########################################################################
    ## end et uStart control parameters
    ###########################################################################

    ###########################################################################
    ## create coxph arguments
    ###########################################################################
    if (hazRestrict=="none") {
        coxX <- x
        coxY <- survival::Surv(y, d)
        coxInt <- NULL  # no intercept columns
    } else if (hazRestrict=="proportional") {
        coxX <- makeBlockDiag(x, p, k)  # Make x into a nk x pk block diagonal matrix
        coxY <- survival::Surv(rep(y, k), rep(d, k))
        ## intercept columns:
        coxInt <- sapply(1:k, function(i) c(rep(0, n*(i-1)), rep(1, n), rep(0, n*(k-i))))
    }
    else {
        coxX <- makeBlockDiag(x, p, k)  # Make x into a nk x pk block diagonal matrix
        coxY <- survival::Surv(rep(y, k), rep(d, k))
        coxInt <- NULL  # no intercept columns
    }
    ###########################################################################
    ## end create coxph arguments
    ###########################################################################

    ###########################################################################
    ## randomly generate uStart (if not passed as argument)
    ###########################################################################
    if (is.null(uStart)) {
        if (uStartMethod=="random") {
            uStart <- getUStartRandom(n, k, uRandomControl$nStart, uRandomControl$p)
        } else if (uStartMethod=="subset") {
            uStart <- getUStartSubset(x, y, d, k, uThresh=0,
                                      npp=uSubsetControl$npp,
                                      nStart=uSubsetControl$nStart,
                                      maxiter=uSubsetControl$maxiter)
        }
    }

    lapply(uStart, function(u) {
        if (nrow(u)!=n || ncol(u)!=k) stop("Wrong dimensions in uStart")
    })
    ###########################################################################
    ## end randomly generate uStart (if not passed as argument)
    ###########################################################################

    ###########################################################################
    ## get list of fitted models (one for each set of uStart values)
    ###########################################################################
    cacLS <- lapply(uStart, function(u) {
        tryCatch({
            cacSingleStart(uStart=u,
                           n=n, p=p, coxX=coxX, coxY=coxY, coxInt=coxInt,
                           uniqY=uniqY, yCt=yCt,
                           x=x, y=y, d=d, k=k,
                           hazRestrict=hazRestrict,
                           maxiter=maxiter, eps=eps, uThresh=0, trace=trace,
                           fixedBeta=fixedBeta, fixedPi=fixedPi)
        }, error = function(e) NULL)
    })

    ###########################################################################
    ## end get list of fitted models (one for each set of uStart values)
    ###########################################################################

    ###########################################################################
    ## return "cacFit" object, including fit with best observed semiparametric
    ## log-likelihood
    ###########################################################################
    loglikLocal <- sapply(cacLS, function(x) if(is.null(x)) NA else x$loglik)
    bestLocalIndex <- which.max(loglikLocal)
    arglist <- list(x=x, y=y, d=d, k=k, hazRestrict=hazRestrict, uStart=uStart,
                    maxiter=maxiter, eps=eps,
                    fixedBeta=fixedBeta, fixedPi=fixedPi)
    cacFit <- c(cacLS[[bestLocalIndex]],
                list(arglist=arglist,
                     bestLocalIndex=bestLocalIndex,
                     loglikLocal=loglikLocal))
    if (keepLocal) cacFit <- c(cacFit, list(localFits=cacLS))
    class(cacFit) <- "cacFit"
    cacFit
}
