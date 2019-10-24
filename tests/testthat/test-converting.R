data("methrix_data")

m2 <- convert_methrix(methrix_data)

test_that("Expected results", {
  #expect_equal(convert_HDF5_methrix(convert_methrix(methrix_data)), methrix_data) #Assays differ in classes (matrix v/s DelayedMatrix)
  expect_equivalent(convert_methrix(convert_HDF5_methrix(m2)), m2)
  })

test_that("Expected errors", {
  expect_error(convert_methrix("not methrix"))
  expect_error(convert_methrix(m2))
  expect_error(convert_HDF5_methrix("not methrix"))
  expect_error(convert_HDF5_methrix(methrix_data))
})
