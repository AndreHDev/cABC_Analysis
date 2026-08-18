#' cABC_plotGG
#' 
#' ggplot2 version matching base R cABC_plot
#'
#' @importFrom grDevices colors
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
#' 
#' @param y_vals Numeric vector of y coordinates of original data points.
#' 
#' @param LineWidth Numeric. Line width for the ABC curve. Default is 3.
#' 
#' @param ShowUniform Logical. If TRUE (default), the uniform reference curve is
#'   drawn in addition to the identity and ABC curves.
#'   
#' @param Plot_title Character string. Title of the plot. Default is "ABC plot".
#' 
#' @details
#' The plot always uses a square coordinate system with both axes ranging from 0 to 1.
#' The diagonal y = 1 - x (equilibrium line) and the identity line y = x
#' are drawn as references. ABC set boundaries (A|B and B|C) are visualized with
#' stars and orthogonal boundary lines. 
#' Shows individual points if they are less then 20.
#' 
#' @return ggplot2 object
#' @keywords internal
cABC_plotGG <- function(CurveData, CleanData, Boundaries, Set_counts, x_vals, y_vals,
                        plot_args = cABC_resolve_plot_args(NULL)) {
  
  # Fall back to full defaults if caller passes NULL or a partial list directly
  plot_args <- cABC_resolve_plot_args(plot_args)
  
  LineWidth         <- plot_args$LineWidth
  ShowUniform       <- plot_args$ShowUniform
  ShowBoundary      <- plot_args$ShowBoundary
  Plot_title        <- plot_args$Plot_title
  UniformColor      <- plot_args$UniformColor
  IdentityColor     <- plot_args$IdentityColor
  EquilibriumColor  <- plot_args$EquilibriumColor
  CurveColor        <- plot_args$CurveColor
  PointColor        <- plot_args$PointColor
  ABoundaryColor    <- plot_args$ABoundaryColor
  BBoundaryColor    <- plot_args$BBoundaryColor
  CBoundaryColor    <- plot_args$CBoundaryColor
  BoundaryLineColor <- plot_args$BoundaryLineColor
  LabelColor        <- plot_args$LabelColor
  LegendTextSize    <- plot_args$LegendTextSize
  Theme             <- plot_args$Theme
  
  Effort <- CurveData$Curve[, 'Effort']
  Yield <- CurveData$Curve[, 'Yield']
  cleaned_data <- CleanData
  
  # Calculate uniform curve
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
  
  # Create base plot
  p <- ggplot2::ggplot() +
    ggplot2::scale_x_continuous(expand = c(0, 0), limits = c(0, 1), 
                                name = 'fraction of data') +
    ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(0, 1),
                                name = 'fraction of sum of largest data') +
    ggplot2::coord_cartesian(
      xlim = c(0, 1),
      ylim = c(0, 1),
      expand = FALSE
    ) +
    ggplot2::coord_fixed(ratio = 1) +
    ggplot2::labs(title = Plot_title)
  
  # Add uniform curve if requested
  if(ShowUniform) {
    uniform_df <- data.frame(x = p_unif, y = ABC_uniform)
    # Clip values to 0 and 1
    uniform_df$x <- pmin(pmax(uniform_df$x, 0), 1)
    uniform_df$y <- pmin(pmax(uniform_df$y, 0), 1)
    p <- p + ggplot2::geom_line(data = uniform_df, ggplot2::aes(x = x, y = y), 
                                color = UniformColor, linewidth = 1)
  }
  
  # Identity line
  p <- p + ggplot2::geom_line(data = data.frame(x = c(0, 1), y = c(0, 1)),
                              ggplot2::aes(x = x, y = y), color = IdentityColor, linewidth = 0.1)
  
  # Equilibrium diagonal
  p <- p + ggplot2::geom_line(data = data.frame(x = c(0, 1), y = c(1, 0)),
                              ggplot2::aes(x = x, y = y), linetype = 'dashed', 
                              color = EquilibriumColor, linewidth = 1)
  
  # ABC curve
  abc_df <- data.frame(x = Effort, y = Yield)
  p <- p + ggplot2::geom_line(data = abc_df, ggplot2::aes(x = x, y = y), 
                              color = CurveColor, linewidth = LineWidth)
  
  # Show data points if less than 20
  if(length(x_vals) < 20) {
    points_df <- data.frame(x = x_vals, y = y_vals)
    p <- p + ggplot2::geom_point(data = points_df, ggplot2::aes(x = x, y = y),
                                 shape = 1, size = 3, color = PointColor, stroke = 1.5)
  }
  
  # Boundary stars
  if(ShowBoundary){
    p <- p + 
      ggplot2::geom_point(data = data.frame(x = Boundaries$A[1], y = Boundaries$A[2]),
                          ggplot2::aes(x = x, y = y), shape = 8, size = 3, color = ABoundaryColor, stroke = 1.5) +
      ggplot2::geom_point(data = data.frame(x = Boundaries$B[1], y = Boundaries$B[2]),
                          ggplot2::aes(x = x, y = y), shape = 8, size = 3, color = BBoundaryColor, stroke = 1.5) +
      ggplot2::geom_point(data = data.frame(x = Boundaries$C[1], y = Boundaries$C[2]),
                          ggplot2::aes(x = x, y = y), shape = 8, size = 3, color = CBoundaryColor, stroke = 1.5)
    
    # Boundary lines
    boundary_segments <- data.frame(
      x = c(0, Boundaries$A[1], 0, Boundaries$C[1]),
      y = c(Boundaries$A[2], 0, Boundaries$C[2], 0),
      xend = c(Boundaries$A[1], Boundaries$A[1], Boundaries$C[1], Boundaries$C[1]),
      yend = c(Boundaries$A[2], Boundaries$A[2], Boundaries$C[2], Boundaries$C[2])
    )
    p <- p + ggplot2::geom_segment(data = boundary_segments, 
                                   ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
                                   color = BoundaryLineColor, linewidth = 1)
    
    # A|B and B|C labels
    if(abs(Boundaries$A[1] - Boundaries$C[1]) > 0.1) {
      p <- p + 
        ggplot2::annotate("text", x = Boundaries$A[1], y = Boundaries$A[2], 
                          label = 'A|B', color = LabelColor, size = 3.5, hjust = -0.2, vjust = 1.5) +
        ggplot2::annotate("text", x = Boundaries$C[1], y = Boundaries$C[2], 
                          label = 'B|C', color = LabelColor, size = 3.5, hjust = -0.2, vjust = 1.5)
    } else {
      p <- p + 
        ggplot2::annotate("text", x = Boundaries$A[1] - 0.05, y = Boundaries$A[2] - 0.03, 
                          label = 'A|B', color = LabelColor, size = 3.5) +
        ggplot2::annotate("text", x = Boundaries$C[1] + 0.025, y = Boundaries$C[2] - 0.025, 
                          label = 'B|C', color = LabelColor, size = 3.5)
    }
    
    # A, B, C set labels with counts
    y_label <- Boundaries$A[2] / 4
    y_count <- y_label - 0.05
    
    # Calculate positions with boundary checks to prevent out-of-bounds labels
    pos_A <- Boundaries$A[1] / 2
    pos_B <- (Boundaries$A[1] + Boundaries$C[1]) / 2
    pos_C <- min((Boundaries$A[1] + Boundaries$C[1]) / 2 + max(abs(Boundaries$A[1] - Boundaries$C[1]), 0.1), 0.95)
    pos_C_count <- min((Boundaries$A[1] + Boundaries$C[1]) / 2 + max(abs(Boundaries$A[1] - Boundaries$C[1]), 0.1) + 0.02, 0.97)
    
    p <- p + 
      # Set A
      ggplot2::annotate("text", x = pos_A, y = y_label,
                        label = 'A', color = LabelColor, size = 9, fontface = 'bold') +
      ggplot2::annotate("text", x = pos_A, y = y_count,
                        label = paste0('n=', Set_counts$nA), color = 'black', size = 2.8) +
      # Set B
      ggplot2::annotate("text", x = pos_B, y = y_label,
                        label = 'B', color = LabelColor, size = 7.5, fontface = 'bold') +
      ggplot2::annotate("text", x = pos_B, y = y_count,
                        label = paste0('n=', Set_counts$nB), color = 'black', size = 2.8) +
      # Set C
      ggplot2::annotate("text", x = pos_C, y = y_label, 
                        label = 'C', color = LabelColor, size = 6.5, fontface = 'bold') +
      ggplot2::annotate("text", x = pos_C_count, y = y_count, 
                        label = paste0('n=', Set_counts$nC), color = 'black', size = 2.8)
  }

  
  # Legend position
  legend_y_start <- if(((Boundaries$A[1] + Boundaries$C[1]) / 2 + 
                        max(abs(Boundaries$A[1] - Boundaries$C[1]), 0.1) + 0.02) < 0.8) {
    0.3
  } else {
    0.5
  }
  
  # Legend annotations
  legend_spacing <- 0.05
  legend_x <- 0.80
  
  legend_entries <- list(
    list(show = ShowBoundary, label = 'set limits', color = BoundaryLineColor),
    list(show = TRUE,         label = 'data',       color = CurveColor),
    list(show = ShowUniform,  label = 'uniform',    color = UniformColor),
    list(show = TRUE,         label = 'identity',   color = IdentityColor)
  )
  legend_entries <- Filter(function(e) e$show, legend_entries)
  
  for (i in seq_along(legend_entries)) {
    e <- legend_entries[[i]]
    p <- p + ggplot2::annotate("text", x = legend_x, y = legend_y_start - (i - 1) * legend_spacing,
                               label = e$label, color = e$color, size = LegendTextSize,
                               hjust = 0, fontface = 'italic')
  }
  
  # Theme
  p <- p + Theme 
  
  return(p)
}