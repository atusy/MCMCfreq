
#' Frequency graph
#'
#' Drawing frequency graph (density plot) and analyzing the out put of
#' \code{\link{cmpt_auto_freq}} or \code{\link{cmpt_freq}}.
#'
#' @export
#' @param data list data of the result of functions; \code{auto_cmpt_freq()} or
#'   \code{cmpt_freq().}
#' @param x_min numeric; the smallest end point of x axis.
#' @param x_max numeric; the largest end point of x axis.
#' @param y_max numeric; the largest end point of y axis.
#' @param uncertainty logical; if \code{TRUE} (default), drawing 95 percent
#'   confident interval of the estimate line.
#' @param rug_plot apply data points' shape of rug plot.
#' \itemize{
#'  \item{0: no plot}
#'  \item{1: thin line plot (default)}
#'  \item{2: thin square plot. This plot should be used for discrete data,
#'  such as integer data.}
#'  \item{3: thick square plot. Also used for discrete data}
#' }
#' @param graphic_theme one of either:
#' \itemize{
#'  \item{0: default theme, with y
#'   axis label, graph frame, and grid.}
#'  \item{1: simple theme, only with x axis
#'   in plot area.}
#' }
#' @param grid one of either (only used if \code{graphic_theme = 1}):
#'  \itemize{
#'   \item{0: no grid}
#'   \item{1: dotted grid (default)}
#'  }
#' @param xlab a text style title for the x axis. Defaults to "Age (Ma)".
#' @param ylab a text style title for the y axis, only used if
#'   \code{graphic_theme = 1}.
#' @param grid_color the color to be used for grid, only used if
#'   \code{graphic_theme = 1}. Defaults to "grey30".
#' @param area_color the color to be used for 95 percent confident interval
#'   area, only used if \code{uncertainty = 1}. Defaults to "grey80".
#' @param line_color the color to be used for estimate line. Defaults to
#'   "black".
#' @param rug_plot_color the color to be used for rug plot. Defaults to rgb(0,
#'   0, 0, alpha = .6).
#' @param hist logical; if \code{TRUE}, drawing histogram in plot area. Defaults
#'   to \code{FALSE}.
#' @param hist_bin integer, larger than 0. Number of histogram bins. Only used
#'   if \code{hist = TRUE}.
#' @param hist_height numeric, between 0 and 1. If \code{hist_height = 1}, then
#'   the height of largest bin is that of plot area. Generally, the height is
#'   \code{hist_height} times. Only used if \code{hist = TRUE}.
#' @param hist_color the color to be used for histgrams. Defaults to rgb(1,
#'   0.62, 0, alpha = .6). Only used if \code{hist = TRUE}.
#' @param analyze logical; if \code{TRUE}, drawing the peak x-coordinate,
#'   largest and smallest data, and number of data. Defaults to \code{FALSE}.
#' @return data.frame of the x-coordinats and height of peaks. Only
#'   returned if \code{analyze =TRUE}.
#'


freq_graph <- function(data,
                       x_min = NA,
                       x_max = NA,
                       y_max = NA,
                       uncertainty = TRUE,
                       rug_plot = 1,
                       graphic_theme = 1,
                       grid = 1,
                       xlab = "Age (Ma)",
                       ylab = "frequency",
                       grid_color = "grey30",
                       area_color = "grey80",
                       line_color = "black",
                       rug_plot_color = rgb(0, 0, 0, alpha = .6),
                       hist = FALSE,
                       hist_bin = as.integer(50),
                       hist_height = 0.2,
                       hist_color = rgb(1, 0.62, 0, alpha = .6),
                       analyze = FALSE){

  #### attentions!  ####
  #-----------------------------------------------------------------------------

  if(is.vector(data) == FALSE) {
    stop("'data' must be vector")
  }
  if(is.logical(uncertainty) == FALSE) {
    stop("'uncertainty' must be logical")
    }
  if(rug_plot != 0 & rug_plot != 1 & rug_plot != 2 & rug_plot !=3 ) {
    stop("'rug_plot' must be 0, 1, 2 or 3")
    }
  if(grid != 1 & grid != 2) {
    stop("'grid' must be 0 or 1")}
  if(graphic_theme != 1 & graphic_theme != 2) {
    stop("'graphic_theme' must be 1 or 2")
    }
  if(is.character(xlab) == F){
    stop("'xlab' must be character")
  }
  if(is.character(ylab) == F){
    stop("'ylab' must be character")
  }
  if(is.character(grid_color) == F){
    stop("'grid_color' must be character")
  }
  if(is.character(area_color) == F){
    stop("'area_color' must be character")
  }
  if(is.character(line_color) == F){
    stop("'line_color' must be character")
  }
  if(is.logical(hist) == F){
    stop("'hist' must be logical")
  }
  if(is.integer(hist_bin) == F | hist_bin < 1){
    stop("'hist_bin' must be integer, > 0")
  }
  if(is.character(hist_color) == F){
    stop("'hist_color' must be character")
  }
  if(0 > hist_height | hist_height >1 | length(hist_height) != 1){
    stop("'hist_height' must be 0 to 1 numeric ")
  }
  if(is.logical(analyze) == F){
    stop("'analyze' must be logical")
  }


  # preperation
  #-----------------------------------------------------------------------------
  p <- data$p
  d <- data$data

  if(is.na(x_min) == T){
    x_min <- min(p$age_X)
  }

  if(is.na(x_max) == T){
    x_max <- max(p$age_X)
  }

  if(x_max <= x_min) {stop("'x_max' must be larger than 'x_min'")}

  if(graphic_theme == 2) {ylab <- ""}


  plot_data <- data.frame(age_X = p$age_X, mean = p$p_mean, p_0025 = p$p_0025, p_0975 = p$p_0975)

  # plot area
  #-----------------------------------------------------------------------------
  if(is.na(y_max) == T){
    y_max <- max(plot_data$p_0975)
    y_max_ceiling <- ceiling(10*y_max)/10
  } else {
    y_max_ceiling <- y_max
  }
  if(y_max < 0) {
    stop("'y_max' must be numeric > 0")
  }

  x <- seq(x_min, x_max, by=(x_max-x_min))
  y <- seq(0, y_max_ceiling, by=y_max_ceiling)


  if(rug_plot == 1 | rug_plot == 2 | rug_plot == 3){
    plot(x, y, ylab=ylab, xlab=xlab, type='n',xaxs="i", axes=FALSE)
  }

  if(rug_plot == 0){
    plot(x, y, ylab=ylab, xlab=xlab, type='n',xaxs="i",yaxs="i", axes=FALSE)
  }

  # rug plot
  #-----------------------------------------------------------------------------

  if(rug_plot == 1){
    for(i in 1:length(d)){
      segments(d[i],0,d[i],-0.2, col=rug_plot_color)
    }
  }

  if(rug_plot == 2){
    width <- max(d) - min(d)
    plot_width <- 1/100 * width

    for(i in 1:length(d)){
      polygon( c(d[i] - plot_width, d[i] + plot_width, d[i] + plot_width,
                 d[i] - plot_width, d[i] - plot_width),
               c(0, 0, -0.1, -0.1, NA),
               col = rug_plot_color,
               lty="blank")
    }
  }

  if(rug_plot == 3){
    width <- max(d) - min(d)
    plot_width <- 1/50 * width

    for(i in 1:length(d)){
      polygon( c(d[i] - plot_width, d[i] + plot_width, d[i] + plot_width,
                 d[i] - plot_width, d[i] - plot_width),
               c(0, 0, -0.1, -0.1, NA),
               col = rug_plot_color,
               lty = "blank")
    }
  }



  abline(a = 0, b = 0)


  # graphic theme
  #-----------------------------------------------------------------------------

  if(graphic_theme == 1){
    if(grid==1){
      axis(side=1,tck=1.0,lty="dotted",col=grid_color, lwd = 0.7)
      axis(side=2,tck=1.0,lty="dotted",col=grid_color, lwd = 0.7)
    }
    box()
  }

  # x axis
  #-----------------------------------------------------------------------------
  #large scale
  axis(side=1, labels=T,tck=-0.025)

  #minimun scale
  xgap <- par("xaxp")
  delt <- (xgap[2] - xgap[1]) / xgap[3] /10
  a <- (xgap[1] - x_min) %/% delt
  b <- (x_max - xgap[2]) %/% delt

  axis(side = 1, labels = F, tck = -0.01,
       xaxp = c(xgap[1] - delt * a, xgap[2] + delt * b, 10 * xgap[3] + a + b))

  # histogram
  #-----------------------------------------------------------------------------

  if(hist == 1){
    bin_width = (x_max - x_min) / hist_bin
    y <- numeric(hist_bin)

    for(j in 1:hist_bin){
      for(i in 1:length(d)){
        if(x_min + (j-1) * bin_width <= d[i] & d[i] < x_min + j * bin_width){
          y[j] = y[j]+1
        }
      }
    }

    y[hist_bin] <- y[hist_bin] + 1

    h <- y_max_ceiling / max(y) * hist_height

    for(j in 1:hist_bin){
      polygon(c(x_min + (j-1) * bin_width, x_min + j * bin_width,
                x_min + j * bin_width, x_min + (j-1) * bin_width,
                x_min + (j-1) * bin_width),
              c(h * y[j], h * y[j], 0, 0, h * y[j]),
              col = hist_color,
              lty = "blank")
    }
  }

  # uncertainty area
  #-----------------------------------------------------------------------------
  if(uncertainty == 1){
    polygon(
      x = c(plot_data$age_X , rev(plot_data$age_X ), plot_data$age_X[1]),
      y = c(plot_data$p_0025, rev(plot_data$p_0975), NA),
      col = rgb(0, 0, 0,alpha = .3),
      lty = "blank"
    )
  }

  # density estimate
  #-----------------------------------------------------------------------------

  lines(plot_data$age_X, plot_data$mean, type = "l", lwd = 2, col = line_color)

  # data analyze
  #-----------------------------------------------------------------------------

  if(analyze == 1){
    mean <- plot_data$mean
    peak <- data.frame(age_X = 0, mean = 0)
    for(i in 2:(length(mean) - 1)){
      if(mean[i-1] < mean[i] & mean[i] > mean[i+1]){
        peak <- rbind(peak, c(plot_data$age_X[i], mean[i]))
      }
    }

    peak <- peak[-1,]
    rownames(peak) <- c()

    points(peak$age_X, peak$mean)

    for(i in 1:length(peak$age_X)){
      age <- signif(peak$age_X[i], digits = 3)
      text = as.character(age)
      text(peak$age_X[i], peak$mean[i] + y_max_ceiling * 0.05, text)
    }

    text2 <- paste("N=", as.character(length(d)), "\n",
                   "smallest_x=", as.character(signif(min(d), digits = 4)), "\n",
                   "largest_x=" , as.character(signif(max(d), digits = 4)))
    text(x_max - (x_max - x_min) * 1/5, 0.9*y_max_ceiling, text2)

    return(peak)
  }
}

