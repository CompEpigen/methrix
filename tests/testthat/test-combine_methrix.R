
library(methrix)
data("methrix_data")

test_that("Combinations, non_hdf5", {
  expect_equivalent(combine_methrix(methrix_data[,1:2], methrix_data[,3:4], by = "col"), methrix_data)
  expect_equivalent(combine_methrix( methrix_data[1:100,], methrix_data[-(1:100),],by = "row"), methrix_data)
})

methrix_data_h5 <- convert_methrix(methrix_data)
test_that("Combinations, hdf5", {
  expect_equivalent(combine_methrix(methrix_data_h5[,1:2], methrix_data_h5[,3:4], by = "col"), methrix_data_h5)
  expect_equivalent(combine_methrix( methrix_data_h5[1:100,], methrix_data_h5[-(1:100),],by = "row"), methrix_data_h5)
})


test_that("Expected errors", {
  expect_error(combine_methrix( methrix_data[1:100,3:4], methrix_data[101:200,1:2], by = "col"), "You have to have the same regions in your datasets.")
  expect_error(combine_methrix( methrix_data[1:100,3:4], methrix_data[50:200,1:2], by = "col"), "You have to have the same regions in your datasets.")
  expect_error(combine_methrix( methrix_data_h5[1:100,3:4], methrix_data_h5[101:200,1:2], by = "col"), "You have to have the same regions in your datasets.")
  expect_error(combine_methrix( methrix_data_h5[1:100,3:4], methrix_data_h5[50:200,1:2], by = "col"), "You have to have the same regions in your datasets.")
  expect_error(combine_methrix( methrix_data, methrix_data, by = "col"), "You have the same samples in your datasets. You need different samples for this merging.")
  expect_error(combine_methrix( methrix_data, methrix_data, by = "row"), "There are overlapping regions in your datasets. This function only takes distinct objects.")
  expect_error(combine_methrix( methrix_data[1:100, 1:3], methrix_data[-(1:100),], by = "row"), "You have different samples in your dataset. You need the same samples in your datasets.")
})
