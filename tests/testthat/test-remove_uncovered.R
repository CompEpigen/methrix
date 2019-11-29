data('methrix_data')


m <- methrix_data
assays(m)$beta[1,] <- NA
assays(m)$cov[1,] <- NA
m2 <- convert_methrix(m)

test_that("Expected results", {
  #expect_failure(remove_uncovered(methrix_data), methrix_data) #This should be a test case for not equal (743 rows v/s 742 rows)
  expect_equivalent(remove_uncovered(m), methrix_data[-1,])
  expect_equivalent(remove_uncovered(m2), convert_methrix(methrix_data[-1,]))
})

test_that("Expected errors", {
  expect_error(remove_uncovered("not methrix"))
})
