#' The 'MCMCfreq' package.
#'
#' @description Package to draw density curve and computation of frequency by MCMC. (Beta Version)
#'
#' @docType package
#' @name MCMCfreq-package
#' @aliases MCMCfreq
#' @useDynLib MCMCfreq, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @import rstantools
#' @importFrom rstan sampling
#' @author Tan Furukawa <rpackagetan@gmail.com>
#' @references
#' \itemize{
#' \item{
#' Stan Development Team (2018). RStan: the R interface to Stan. R package version 2.18.2. http://mc-stan.org
#' }
#' \item{
#' Tsujimori, T. (1998). Geology of the Osayama serpentinite melange in the central Chugoku Mountains,
#'  southwestern Japan: 320 Ma blueschist-bearing serpentinite melange beneath the Oyama ophiolite.
#'  Jour. Geol. Soc. Japan, 104, 213-231.
#' }
#' \item{
#' Riihimäki, J., & Vehtari, A. (2014). Laplace approximation for logistic Gaussian process density
#' estimation and regression. Bayesian analysis, 9(2), 425-448.
#' }
#' }
NULL
