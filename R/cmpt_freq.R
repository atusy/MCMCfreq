
#' Bayesian Gaussian Process regression with Stan: computation of frequency by MCMC
#'
#' The frequency is computed by Gaussian Process (GP) regression. In this
#' computation, only one-dimensional data can be analyzed, but accurate
#' analysis is possible. The frequency is expressed by parameter 'p' distrubution,
#' and p is defined by p = softmax(f), f ~ multi_normal(0, K).
#' K is computed by Gaussian Kernel; K(x, x') = sigma^2 * exp((x - x')^2 / rho^2).
#' For explain briefly, sigma and rho govern the smoothness properities.
#' Sigma is related to amplitude of frequency plot, and rho is how fast
#' the correlation decrease (it means that if rho become smaller, the possibility is
#' increasing for having a completely different frequency in distant two points of x).
#' sigma and rho should be determined by the user, however, there is no clue in
#'  many cases. Then, function \code{cmpt_freq()} give us a hint; WAIC.
#' WAIC is one of an estimator of the relative quality of statistical models
#' for a given set of data. If WAIC is smaller, it can be said that it is a
#' better model.
#'
#' @export
#' @param data a vector of one dimentional data.
#' @param rho numeric >0; how fast the correlation decrease.
#' (See \strong{Description})
#' @param sigma numeric >0; amplitude of frequency plot.
#' (See \strong{Description})
#' @param delta numeric beetween 0 and 1 (defaults to 1/25); the size of
#' no-data area in frequency plot.
#'  In order to stabilize computation, both ends of the graph have no data
#' points. The size is expressed by \code{delta} times size of plot area
#'  in the x-axis direction.
#' @param M integer; degree of discretization computation. x axis of
#' frequency plot is divided into M grids.
#'  Defaults to \code{M = 50}, and larger M give accurate frequency plot,
#'   but the computation takes much time.
#' @param stan_seed see \code{?stan} in \code{library(rstan)}
#' @param stan_chains see \code{?stan} in \code{library(rstan)}
#' @param stan_warmup see \code{?stan} in \code{library(rstan)}
#' @param stan_thin see \code{?stan} in \code{library(rstan)}
#' @param stan_iter see \code{?stan} in \code{library(rstan)}
#' @param stan_max_treedepth see \code{?stan} in \code{library(rstan)}
#' @return list data having follow 5 components:
#' \itemize{
#'  \item{\code{$WAIC}: See \strong{Description}}
#'  \item{\code{$p}: data.frame of the MCMC result. \code{$age_X} is x
#'  coodinates and \code{$p_mean} is estimate. \code{$p_0025} means that
#'  the probability that ture frequency is smaller than \code{$p_0025}
#'  is 0.025. \code{$p_025}, \code{$p_075} and \code{$p_0975} is defined
#'  in the same way as \code{$p_0025}. \code{$p_n_eff} and \code{$p_Rhat} is
#'  index to check convergence. In this function, if \code{$p_0025 < ##########} and \code{$p_Rhat < ########}, then consider the MCMC to converge}
#'  \item{\code{$data}: same as the input data, used in \code{\link{freq_graph}}}
#'  \item{\code{$rho}: same as the input rho}
#'  \item{\code{$sigma}: same as the input sigma}
#' }

cmpt_freq <- function(data,
                      rho,
                      sigma,
                      delta = 1/25,
                      M = as.integer(50),
                      stan_seed = as.integer(1234),
                      stan_chains = as.integer(3),
                      stan_warmup = as.integer(300),
                      stan_thin = as.integer(1),
                      stan_iter = as.integer(1300),
                      stan_max_treedepth = as.integer(10)){
  ##### attention ! #####
  #-----------------------------------------------------------------------------
  if(is.vector(data) == FALSE) {
    stop("'data' must be vector")
  }
  if(is.integer(M) == FALSE | M < 0) {
    stop("'M' must be integer")
  }
  if(is.numeric(delta) == FALSE | delta < 0 | delta > 1) {
    stop("'delta' must be numeric between 0 and 1")
  }
  if(is.numeric(rho) == FALSE | delta < 0) {
    stop("'rho' must be numeric > 0")
  }
  if(is.numeric(sigma) == FALSE | delta < 0) {
    stop("'sigma' must be numeric > 0")
  }
  if(is.integer(stan_seed) == FALSE | stan_seed < 0) {
    stop("'stan_seed' must be integer")
  }
  if(is.integer(stan_chains) == FALSE | stan_chains < 0) {
    stop("'stan_chains' must be integer")
  }
  if(is.integer(stan_warmup) == FALSE | stan_warmup < 0) {
    stop("'stan_warmup' must be integer")
  }
  if(is.integer(stan_thin) == FALSE | stan_thin < 0) {
    stop("'stan_thin' must be integer")
  }
  if(is.integer(stan_iter) == FALSE | stan_iter < 0) {
    stop("'stan_iter' must be integer")
  }
  if(is.integer(stan_max_treedepth) == FALSE | stan_max_treedepth < 0) {
    stop("'stan_max_treedepth' must be integer")
  }

  #package etc
  #-----------------------------------------------------------------------------

  auto_write=T
  options(mc.cores=parallel::detectCores())

  #Data preparation
  #-----------------------------------------------------------------------------
  d <- data
  N <- length(d)
  max_d <- max(d)
  min_d <- min(d)
  #the definiton of various

  #     |<---------------  1 --------------->|
  #     delta                            delta            normalized_age
  #     |   |<------ 1 - 2*delta ------->|   |
  #------------------------------------------------------> age
  #     |   | <= min_d          max_d => |   |
  #     Delta                            Delta            real_age
  #     | <= Min_d                  Max_d => |

  Delta <- (max_d-min_d)/(1-2*delta) * delta
  Min_d <- min_d - Delta
  Max_d <- max_d + Delta

  #nd is a normalized value of d
  nd <-  (d - Min_d)/(Max_d-Min_d)
  #X
  X <- numeric(M)
  for(i in 1:M){X[i] <- 1/(2*M)+(i-1)/M}
  #y
  y <- numeric(M)
  for(j in 1:M){
    for(i in 1:N){
      if((j-1)*(1/M)<=nd[i] & nd[i]<=j*(1/M)){
        y[j] = y[j]+1
      }
    }
  }


  #function for finding WAIC
  #-----------------------------------------------------------------------------

  dlist <- list(M = M, x = X*100, y = y, rho=rho^2, sigma=sigma^2)
  fit <- rstan::sampling(stanmodels$density,
                  data = dlist,
                  seed = stan_seed,
                  chains = stan_chains,
                  warmup = stan_warmup,
                  iter = stan_iter,
                  thin = stan_thin,
                  control = list(max_treedepth = stan_max_treedepth))
  print("computing p value.....")
  p <- data.frame(p_mean  = numeric(M),
                  p_0025  = numeric(M),
                  p_025   = numeric(M),
                  p_075   = numeric(M),
                  p_0975  = numeric(M),
                  p_n_eff = numeric(M),
                  p_Rhat  = numeric(M))

  MCMC_header <- c("mean" ,
                   "2.5%" ,
                   "25%"  ,
                   "75%"  ,
                   "97.5%",
                   "n_eff",
                   "Rhat" )

  for(j in 1:length(MCMC_header)){
    for(i in 1:M){
      name <- sprintf("p[%d]",i)
      p[i,j] <- rstan::summary(fit)$summary[name,MCMC_header[j]]
    }
  }

  if(prod(p$p_Rhat < 1.1) == 0 ){
    warning("MCMC is not convergent (Rhat > 1.1)")
  }
  cat("computing WAIC.....")
  ##WAIC
  q <- rstan::extract(fit)$p
  q2 <- colMeans(q)
  q3 <- log(q)
  q4 <- apply(q3, 2, var)

  lppd <- sum(y*log(q2))
  pwaic <- sum(y*q4)
  WAIC <- -2*(lppd-pwaic)

  age_X <- Min_d + (Max_d-Min_d) * X
  p <- cbind(data.frame(age_X=age_X),p)

  answer <- list(WAIC = WAIC, p = p, data = d, rho = rho, sigma = sigma)

  return(answer)

}

