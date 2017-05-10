makeBlockDiag <- function(x, p, k)
{
    bdRow <- function(i) {
        leftBlock <- c(rep(list(0), p*(i-1)))
        midBlock <- list(x)
        rightBlock <- rep(list(0), p*(k-i))
        blockList <- c(leftBlock, midBlock, rightBlock)
        do.call(cbind, blockList)
    }
    blockX <- do.call(rbind, lapply(1:k, bdRow))
    blockX
}
