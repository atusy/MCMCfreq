% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/freq_graph.R
\name{freq_graph}
\alias{freq_graph}
\title{Frequency graph}
\usage{
freq_graph(data, x_min = NA, x_max = NA, y_max = NA,
  uncertainty = TRUE, rug_plot = 1, graphic_theme = 1, grid = 1,
  xlab = "Age (Ma)", ylab = "frequency", grid_color = "grey30",
  area_color = "grey80", line_color = "black",
  rug_plot_color = rgb(0, 0, 0, alpha = 0.6), hist = FALSE,
  hist_bin = as.integer(50), hist_height = 0.2, hist_color = rgb(1,
  0.62, 0, alpha = 0.6), analyze = FALSE)
}
\arguments{
\item{data}{list data of the result of functions; \code{auto_cmpt_freq()} or
\code{cmpt_freq().}}

\item{x_min}{numeric; the smallest end point of x axis.}

\item{x_max}{numeric; the largest end point of x axis.}

\item{y_max}{numeric; the largest end point of y axis.}

\item{uncertainty}{logical; if \code{TRUE} (default), drawing 95 percent
confident interval of the estimate line.}

\item{rug_plot}{apply data points' shape of rug plot.
\itemize{
 \item{0: no plot}
 \item{1: thin line plot (default)}
 \item{2: thin square plot. This plot should be used for discrete data,
 such as integer data.}
 \item{3: thick square plot. Also used for discrete data}
}}

\item{graphic_theme}{one of either:
\itemize{
 \item{1: default theme, with y
  axis label, graph frame, and grid.}
 \item{2: simple theme, only with x axis
  in plot area.}
}}

\item{grid}{one of either (only used if \code{graphic_theme = 1}):
\itemize{
 \item{0: no grid}
 \item{1: dotted grid (default)}
}}

\item{xlab}{a text style title for the x axis. Defaults to "Age (Ma)".}

\item{ylab}{a text style title for the y axis, only used if
\code{graphic_theme = 1}.}

\item{grid_color}{the color to be used for grid, only used if
\code{graphic_theme = 1}. Defaults to "grey30".}

\item{area_color}{the color to be used for 95 percent confident interval
area, only used if \code{uncertainty = 1}. Defaults to "grey80".}

\item{line_color}{the color to be used for estimate line. Defaults to
"black".}

\item{rug_plot_color}{the color to be used for rug plot. Defaults to rgb(0,
0, 0, alpha = .6).}

\item{hist}{logical; if \code{TRUE}, drawing histogram in plot area. Defaults
to \code{FALSE}.}

\item{hist_bin}{integer, larger than 0. Number of histogram bins. Only used
if \code{hist = TRUE}.}

\item{hist_height}{numeric, between 0 and 1. If \code{hist_height = 1}, then
the height of largest bin is that of plot area. Generally, the height is
\code{hist_height} times. Only used if \code{hist = TRUE}.}

\item{hist_color}{the color to be used for histgrams. Defaults to rgb(1,
0.62, 0, alpha = .6). Only used if \code{hist = TRUE}.}

\item{analyze}{logical; if \code{TRUE}, drawing the peak x-coordinate,
largest and smallest data, and number of data. Defaults to \code{FALSE}.}
}
\value{
data.frame of the x-coordinats and height of peaks. Only
  returned if \code{analyze =TRUE}.
}
\description{
Drawing frequency graph (density plot) and analyzing the out put of
\code{\link{cmpt_auto_freq}} or \code{\link{cmpt_freq}}.
}
\examples{
d <- Osayama
e <- cmpt_freq(d,rho=3,sigma=2)
freq_graph(e)
}
