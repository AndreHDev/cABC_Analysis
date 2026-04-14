#' cABC Plot
#'
#' Draws an ABC curve together with identity and optional uniform reference curves,
#' ABC set boundaries (A-B and B-C), labels, and point counts.
#' 
#' @importFrom graphics par plot lines points legend axis box
#' 
#' @param CurveData Data about the ABC Curve as returned by ABC_curve
#' 
#' @param CleanData Clean original input data.
#' 
#' @param Boundaries A list with numeric vectors A, B, and C,
#'   each of length 2, giving the x/y coordinates of the ABC boundaries.
#'
#' @param Set_counts A list with elements nA, nB, nC giving the
#'   number of observations in sets A, B, and C.
#'
#' @param x_vals Numeric vector of x coordinates of original data points.
#' @param y_vals Numeric vector of y coordinates of original data points.
#'
#' @param LineType Integer. If 0 (default), the ABC curve is drawn as a line.
#'
#' @param LineWidth Numeric. Line width for the ABC curve. Default is 3.
#'
#' @param ShowUniform Logical. If TRUE (default), the uniform reference curve is
#'   drawn in addition to the identity and ABC curves.
#'
#' @param Plot_title Character string. Title of the plot. Default is "ABC plot".
#'
#' @param defaultAxes Logical. If TRUE (default FALSE), base R axes are drawn by plot().
#'   If FALSE, custom axes with ticks at 0–1 in steps of 0.1 are drawn.
#'
#' @param ResetPlotDefaults Logical. If TRUE (default), the original par()
#'   settings are restored after plotting.
#'   
#' @details
#' The plot always uses a square coordinate system with both axes ranging from 0 to 1.
#' The diagonal y = 1 - x (equilibrium line) and the identity line y = x
#' are drawn as references. ABC set boundaries (A|B and B|C) are visualized with
#' stars and orthogonal boundary lines.
#' 
#' @return Base R Plot
#' @noRd
cABC_plot <- function(CurveData, CleanData, Boundaries, Set_counts, x_vals, y_vals, LineType=0, LineWidth=3, ShowUniform=TRUE, 
                      Plot_title='ABC plot', defaultAxes = FALSE, ResetPlotDefaults=TRUE) {
  
  Effort <- CurveData$Curve[, 'Effort']
  Yield <- CurveData$Curve[, 'Yield']
  cleaned_data <- CleanData
  
  # Colors and labels
  colors <- c('blue', colors()[452], 'green', colors()[175])
  labels <- c(expression(italic("data")), expression(italic("identity")), 
              expression(italic("uniform")), '')
  
  # Set square plotting region
  par(pty = "s")

  # Uniform curve
  p_unif <- seq(0, 1, by = 0.01)
  if(!is.null(cleaned_data) && length(cleaned_data) > 0) {
    A <- min(cleaned_data, na.rm = TRUE)
    MaxX <- max(cleaned_data, na.rm = TRUE)
    if(A == MaxX) { A <- 0; MaxX <- 1 }
  } else {
    A <- 0; MaxX <- 1
  }
  B <- MaxX - A
  ABC_uniform <- (-0.5 * B * p_unif^2 + MaxX * p_unif) / (A + 0.5 * B)
  
  # Create base plot and plot uniform distribution if wanted
  if(ShowUniform) {
    plot(p_unif, ABC_uniform, type = 'l', col = colors[3], lwd = 1,
         xlim = c(0,1), ylim = c(0,1), xaxs = 'i', yaxs = 'i',
         xlab = 'fraction of data', ylab = 'fraction of sum of largest data',
         main = Plot_title, axes = defaultAxes)
    # Identity distribution
    lines(c(0,1), c(0,1), lty = 1, lwd = 0.1, col = colors()[452])
  } else {
    if(LineType == 0) {
      plot(Effort, Yield, xlim = c(0,1), ylim = c(0,1), xaxs = 'i', yaxs = 'i',
           xlab = 'fraction of data', ylab = 'fraction of sum of largest data',
           type = 'l', lwd = LineWidth, col = colors[1], main = Plot_title, axes = defaultAxes)
    } else {
      plot(Effort, Yield, xlim = c(0,1), ylim = c(0,1), asp = 1, xaxs = 'i', yaxs = 'i',
           pch = LineType, lwd = LineWidth, col = colors[1], main = Plot_title, axes = defaultAxes)
    }
  }
  
  # Equilibrium diagonal
  lines(c(0, 1), c(1, 0), lty = 2, lwd = 1, col = colors()[175])
  
  # ABC curve + data points + stars + lines + annotations 
  if(LineType == 0) {
    lines(Effort, Yield, lwd = LineWidth, col = colors[1])
  } else {
    points(Effort, Yield, pch = LineType, lwd = LineWidth, col = colors[1])
  }
  
  # Show data points if less then 20
  if(length(x_vals) < 20){
    points(x_vals, y_vals, pch = 1, lwd = 1.5, col = 'blue', cex = 1.5)
  }
  
  # Boundary and special points
  points(Boundaries$A[1], Boundaries$A[2], pch=8, lwd=1.5, col='red', cex=1.5)
  points(Boundaries$B[1], Boundaries$B[2], pch=8, lwd=1.5, col='green', cex=1.5)
  points(Boundaries$C[1], Boundaries$C[2], pch=8, lwd=1.5, col='blue', cex=1.5)
  
  # Boundary lines
  lty <- 1
  points(c(0, Boundaries$A[1]), c(Boundaries$A[2], Boundaries$A[2]), type='l', col='red', lty=lty, lwd=1)
  points(c(0, Boundaries$C[1]), c(Boundaries$C[2], Boundaries$C[2]), type='l', col='red', lty=lty, lwd=1)
  points(c(Boundaries$A[1], Boundaries$A[1]), c(0, Boundaries$A[2]), col='red', type='l', lty=lty, lwd=1)
  points(c(Boundaries$C[1], Boundaries$C[1]), c(0, Boundaries$C[2]), col='red', type='l', lty=lty, lwd=1)
  
  # Boundary point annotations
  colors_anno <- c('black', 'red', 'blue', 'green', colors()[452], 'red')
  if(abs(Boundaries$A[1] - Boundaries$C[1]) > 0.1) {
    plotrix::thigmophobe.labels(Boundaries$A[1], Boundaries$A[2], 'A|B', col=colors_anno[2], cex=1)
    plotrix::thigmophobe.labels(Boundaries$C[1], Boundaries$C[2], 'B|C', col=colors_anno[2], cex=1)
  } else {
    plotrix::thigmophobe.labels(Boundaries$A[1]-0.05, Boundaries$A[2], 'A|B', col=colors_anno[2], cex=1)
    plotrix::thigmophobe.labels(Boundaries$C[1]+0.025, Boundaries$C[2]+0.025, 'B|C', col=colors_anno[2], cex=1)
  }
  
  # A,B,C set label
  plotrix::thigmophobe.labels(Boundaries$A[1]/2, Boundaries$A[2]/4, 'A', col=colors_anno[6], cex=2.6)
  plotrix::thigmophobe.labels((Boundaries$A[1]+Boundaries$C[1])/2, Boundaries$A[2]/4, 'B', col=colors_anno[6], cex=2.1)
  plotrix::thigmophobe.labels((Boundaries$A[1]+Boundaries$C[1])/2 + max(abs(Boundaries$A[1]-Boundaries$C[1]), 0.1), 
                              Boundaries$A[2]/4, 'C', col=colors_anno[6], cex=1.8)
  
  # Amount of points per set
  plotrix::thigmophobe.labels(Boundaries$A[1]/2, Boundaries$A[2]/4-0.05, paste0('n=', Set_counts$nA), col='black', cex=0.8)
  plotrix::thigmophobe.labels((Boundaries$A[1]+Boundaries$C[1])/2, Boundaries$A[2]/4-0.05, paste0('n=', Set_counts$nB), col='black', cex=0.8)
  plotrix::thigmophobe.labels((Boundaries$A[1]+Boundaries$C[1])/2 + max(abs(Boundaries$A[1]-Boundaries$C[1]), 0.1) + 0.02, 
                              Boundaries$A[2]/4-0.05, paste0('n=', Set_counts$nC), col='black', cex=0.8)
  
  # Legend + box
  legend_colors <- c('black', 'red', 'blue', 'green', colors()[452], 'red')
  plot_labels <- c('', expression(italic("set limits")), expression(italic("data")), 
                   expression(italic("uniform")), expression(italic("identity")))
  legend_pos <- if(((Boundaries$A[1]+Boundaries$C[1])/2 + max(abs(Boundaries$A[1]-Boundaries$C[1]), 0.1) + 0.02) < 0.8) {
    'bottomright' } else { 'right' }
  legend(legend_pos, legend=plot_labels, text.col=legend_colors, bty="n", y.intersp=0.8)
  
  # Axis
  if(!defaultAxes){
    axis(1, xlim=c(0,1), col="black", las=1, at=seq(0, 1, 0.1))
    axis(2, ylim=c(0,1), col="black", las=1, at=seq(0, 1, 0.1))
  }
  
  # Coordinate system border
  box(col='grey')
  
  # Reset plot if wanted
  if(ResetPlotDefaults) par(par(no.readonly = TRUE))
  
  invisible(list(ABCx = Effort, ABCy = Yield, x_vals = x_vals, y_vals = y_vals))
}
