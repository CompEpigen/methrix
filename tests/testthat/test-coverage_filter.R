
library(methrix)
data("methrix_data")

m1 <- methrix_data
m2 <- convert_methrix(methrix_data)

test_that("Nothing filtered", {
  expect_equivalent(coverage_filter(m1, cov_thr = 1, min_samples = 1), m1)
  expect_equivalent(coverage_filter(m2, cov_thr = 1, min_samples = 1), m2)
})


#Expencted length after filtering.

cov <- get_matrix(m1, type="C")
nbr <- sum(apply(cov, 1, function(x) any(x>=10, na.rm=TRUE)), na.rm=TRUE)

test_that("Expected filtering", {
  expect_equal(nrow(coverage_filter(m1, cov_thr = 10, min_samples = 1)), nbr)
  expect_equal(nrow(coverage_filter(m2, cov_thr = 10, min_samples = 1)), nbr)
  expect_equal(nrow(coverage_filter(m1, cov_thr = 200, min_samples = 1)), 0)
  expect_equal(nrow(coverage_filter(m2, cov_thr = 200, min_samples = 1)), 0)
  expect_equal(nrow(coverage_filter(m1, cov_thr = 10, min_samples = 5)), 0)
  expect_equal(nrow(coverage_filter(m2, cov_thr = 10, min_samples = 5)), 0)
})


test_that("Expected errors", {
  expect_error(coverage_filter("not methrix", cov_thr = 10, min_samples = 5), "A valid methrix object needs to be supplied.")
  expect_error(coverage_filter(m1, cov_thr = "not_number", min_samples = 5), "cov_thr is not numeric.")
  expect_error(coverage_filter(m1, cov_thr = 10, min_samples = "nn"), "min_samples and prop_samples variables are not numeric.")

})

