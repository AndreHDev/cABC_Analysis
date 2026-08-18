library(testthat)
library(cABCanalysis)

# ----------------- Tests for plotArgs / cABC_plot_style -----------------

test_that("NULL plotArgs resolves to defaults", {
  res <- cABC_resolve_plot_args(NULL)
  defaults <- cABC_default_plot_args()
  expect_equal(res, defaults)
})

test_that("partial list is merged onto defaults", {
  res <- cABC_resolve_plot_args(list(LineWidth = 2))
  defaults <- cABC_default_plot_args()
  
  expect_equal(res$LineWidth, 2)
  expect_equal(res$CurveColor, defaults$CurveColor)
  expect_equal(length(res), length(defaults))
})

test_that("unrecognized names trigger a warning and are dropped", {
  expect_warning(
    res <- cABC_resolve_plot_args(list(LineWidth = 2, Foo = "bar")),
    "Foo"
  )
  expect_false("Foo" %in% names(res))
  expect_equal(res$LineWidth, 2)
})

test_that("non-list plotArgs errors", {
  expect_error(
    cABC_resolve_plot_args("not a list"),
    "must be a list"
  )
})

test_that("cABC_plot_style only overrides supplied options", {
  style <- cABC_plot_style(CurveColor = "darkred")
  defaults <- cABC_default_plot_args()
  
  expect_equal(style$CurveColor, "darkred")
  expect_equal(style$LineWidth, defaults$LineWidth)
  expect_equal(style$Theme, defaults$Theme)
})

test_that("cABC_plot_style with no arguments equals defaults", {
  style <- cABC_plot_style()
  defaults <- cABC_default_plot_args()
  expect_equal(style, defaults)
})

test_that("cABC_plot_style accepts a ggplot2 theme", {
  style <- cABC_plot_style(Theme = ggplot2::theme_minimal())
  expect_s3_class(style$Theme, "theme")
})

# ----------------- Integration with cABC_analysis -----------------

test_that("cABC_analysis accepts plotArgs and produces a ggplot object", {
  data <- c(10, 8, 6, 4, 2, 1)
  res <- cABC_analysis(data, PlotIt = TRUE, useGGPlot = TRUE,
                       plotArgs = cABC_plot_style(LineWidth = 2, CurveColor = "darkred"))
  
  expect_s3_class(res$Plot, "ggplot")
})

test_that("cABC_analysis warns on unrecognized plotArgs names", {
  data <- c(10, 8, 6, 4, 2, 1)
  expect_warning(
    res <- cABC_analysis(data, PlotIt = TRUE, useGGPlot = TRUE,
                         plotArgs = list(NotARealOption = 5)),
    "NotARealOption"
  )
  expect_s3_class(res$Plot, "ggplot")
})

test_that("cABC_analysis with PlotIt = FALSE ignores plotArgs entirely", {
  data <- c(10, 8, 6, 4, 2, 1)
  res <- cABC_analysis(data, PlotIt = FALSE, plotArgs = list(LineWidth = 2))
  expect_null(res$Plot)
})