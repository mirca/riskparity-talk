load("Sigma_mu.RData")
library(riskParityPortfolio)
library(quadprog)
library(ggplot2)

max_sharpe_ratio <- function(Sigma, mu) {
    N <- ncol(Sigma)
    if (all(mu <= 1e-8))
        return(rep(0, N))
    Dmat <- 2 * Sigma
    Amat <- diag(N)
    Amat <- cbind(mu, Amat)
    bvec <- c(1, rep(0, N))
    dvec <- rep(0, N)
    res <- solve.QP(Dmat = Dmat, dvec = dvec, Amat = Amat, bvec = bvec, meq = 1)
    w <- res$solution
    return(w/sum(w))
}

rpp_vanilla <- riskParityPortfolio(Sigma)
w_all <- cbind("risk parity" = rpp_vanilla$w, "Markowitz" = max_sharpe_ratio(Sigma, mu))

gr = .5 * (1 + sqrt(5))
setEPS()
postscript("markowitz-rpp-comparison.ps", height = 5, width = 3 * gr)
barplotPortfolioRisk(w_all, Sigma) + theme(legend.position="top")
dev.off()
