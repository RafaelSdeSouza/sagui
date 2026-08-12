test_that("Torch availability requires a loadable runtime", {
  expect_false(sagui:::.torch_runtime_available(
    namespace_available = FALSE,
    runtime_check = function() stop("must not run")
  ))

  expect_false(sagui:::.torch_runtime_available(
    namespace_available = TRUE,
    runtime_check = function() stop("Lantern is not loaded")
  ))

  expect_true(sagui:::.torch_runtime_available(
    namespace_available = TRUE,
    runtime_check = function() TRUE
  ))
})

test_that("distance calculation falls back to base R", {
  x <- matrix(
    c(0, 0, 1, 0, 1, 1),
    ncol = 2,
    byrow = TRUE
  )

  observed <- sagui:::.compute_distance_matrix(x, torch_available = FALSE)
  expected <- stats::dist(x)

  expect_s3_class(observed, "dist")
  expect_equal(as.numeric(observed), as.numeric(expected))
})
