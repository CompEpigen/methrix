

library(methrix)
data("methrix_data")

m1 <- methrix_data
m2 <- convert_methrix(methrix_data)

num_na <- sum(is.na(get_matrix(m1, type = "C")[,1]))
q_99 <- quantile(get_matrix(m1, type = "C")[,1], probs = 0.99, na.rm=TRUE)
num_larger <- sum(get_matrix(m1, type = "C")[,1]>q_99, na.rm = TRUE)


mask_methrix(m1,high_quantile = .99)

test_that("Expected_results", {
  expect_equal(sum(is.na(get_matrix(mask_methrix(m1, low_count = 1, high_quantile = .99), type="C")[,1])), num_na+num_larger)
  expect_equal(sum(is.na(get_matrix(mask_methrix(m2, low_count = 1, high_quantile = .99), type="C")[,1])), num_na+num_larger)
})

test_that("Expected errors", {
  expect_error(mask_methrix("not methrix"), "A valid methrix object needs to be supplied.")
  expect_error(mask_methrix(m1, low_count = "not a number", high_quantile = .99), "low_count must be a numeric value.")
  expect_error(mask_methrix(m1, low_count = 1, high_quantile = "not valid"), "High quantile should be between 0 and 1.")
  expect_equal(nrow(remove_uncovered(mask_methrix(m1, low_count = 300, high_quantile = .99))), 0)
})
