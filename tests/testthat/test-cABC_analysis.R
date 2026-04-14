library(testthat)
library(cABCanalysis)

# ----------------- Tests for cABC_analysis -----------------

test_that("identical values", {
  test1 <- rep(1, 13)
  expect_warning(
    res <- cABC_analysis(test1, PlotIt = FALSE),
    "All data values are identical"
  )
  expect_setequal(unname(res$Aind), 1:13)
})

test_that("character input", {
  test2 <- c("10", 5, "x", 3, "0", 0)
  expect_warning(
    res <- cABC_analysis(test2, PlotIt = FALSE),
    "3 of 6 items are positive"
  )
  expect_setequal(unname(res$Aind), 1)
  expect_setequal(unname(res$Bind), 2)
  expect_setequal(unname(res$Cind), 3:6)
})

test_that("a lot of zeros", {
  test3 <- c(0,10,0,9,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
  warnings <- NULL
  res <- withCallingHandlers(
    cABC_analysis(test3, PlotIt = FALSE),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  
  # Check multiple warnings
  expect_true(any(grepl("3 of 23 items are positive", warnings, fixed = TRUE)))
  expect_true(any(grepl("duplicate value(s) spanning multiple classes", warnings, fixed = TRUE)))
})

test_that("AB duplicate values", {
  test4 <- c(3,1,3,1,1,3)
  expect_warning(
    res <- cABC_analysis(test4, PlotIt = FALSE),
    "duplicate value\\(s\\) spanning multiple classes"
  )
  expect_setequal(unname(res$Aind), c(1,3,6))
  expect_equal(unname(res$Bind), integer(0))
  expect_setequal(unname(res$Cind), c(2,4,5))
})

test_that("single high value, uniform rest", {
  test5 <- c(60, rep(10, 14))
  expect_warning(
    res <- cABC_analysis(test5, PlotIt = FALSE),
    "duplicate value\\(s\\) spanning multiple classes"
  )
  expect_equal(unname(res$Aind), 1)
  expect_equal(unname(res$Bind), integer(0))
  expect_equal(unname(res$Cind), 2:15)
})

test_that("BC duplicate values", {
  test6 <- c(10,10,10,5,5,5,1,1,1)
  expect_warning(
    res <- cABC_analysis(test6, PlotIt = FALSE),
    "duplicate value\\(s\\) spanning multiple classes"
  )
  expect_setequal(unname(res$Aind), c(1,2,3))
  expect_setequal(unname(res$Bind), c(4,5,6))
  expect_setequal(unname(res$Cind), c(7,8,9))
})

test_that("linear data", {
  test7 <- 1:10
  expect_warning(
    res <- cABC_analysis(test7, PlotIt = FALSE),
    NA  # Expect no warnings
  )
})

test_that("high value point", {
  test8 <-  c(100,10,9,8,7,6,5,4,3,2,1)
  res <- cABC_analysis(test8, PlotIt = FALSE)
  expect_equal(unname(res$Aind), 1)
})

test_that("uniform but not identical data", {
  set.seed(123)
  test9 <- 100 + runif(100, min = -0.0002, max = 0.0002)
  expect_warning(
    res <- cABC_analysis(test9, PlotIt = FALSE),
    NA  # Expect no warnings
  )
})

test_that("uniform but not identical data, large scale", {
  set.seed(456)
  test10 <- 100 + runif(10000, min = -0.0002, max = 0.0002)
  expect_warning(
    res <- cABC_analysis(test10, PlotIt = FALSE),
    NA  # Expect no warnings
  )
})

test_that("ABC curve under uniform at the beginning", {
  test11 <- c(10,10,1,1,1)
  expect_warning(
    res <- cABC_analysis(test11, PlotIt = FALSE),
    "duplicate value\\(s\\) spanning multiple classes"
  )
})

test_that("NULL input", {
  expect_error(cABC_analysis(NULL))
})

test_that("single number input", {
  expect_warning(
    res <- cABC_analysis(42, PlotIt = FALSE),
    "Only one data point remains after filtering"
  )
  expect_equal(unname(res$Aind), 1)
  expect_equal(unname(res$Bind), integer(0))
  expect_equal(unname(res$Cind), integer(0))
})

test_that("single number input vector", {
  expect_warning(
    res <- cABC_analysis(c(42), PlotIt = FALSE),
    "Only one data point remains after filtering"
  )
  expect_equal(unname(res$Aind), 1)
  expect_equal(unname(res$Bind), integer(0))
  expect_equal(unname(res$Cind), integer(0))
})

test_that("input with NA values", {
  test12 <- c(10, NA, 5, 7, NA, 3)
  expect_warning(
    res <- cABC_analysis(test12, PlotIt = FALSE),
    "4 of 6 items are positive"
  )
  expect_setequal(unname(res$Aind), c(1,4))
  expect_setequal(unname(res$Bind), 3)
  expect_setequal(unname(res$Cind), c(2,5,6))
})

test_that("<= then 3 values", {
  test13 <- c(10,4)
  expect_warning(
    res <- cABC_analysis(test13, PlotIt = FALSE),
    "Extremly small dataset"
  )
})

test_that("<= then 3 values with zeros", {
  test14 <- c(0,10,0,4,0,0)
  expect_warning(
    res <- cABC_analysis(test14, PlotIt = FALSE),
    "2 of 6 items are positive"
  )
})

test_that("duplicate value spanning all three sets", {
  test15 <- c(1, rep(10, 15))
  expect_warning(
    res <- cABC_analysis(test15, PlotIt = FALSE),
    "duplicate value\\(s\\) spanning multiple classes"
  )
})
