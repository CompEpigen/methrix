data('methrix_data')
m2 <- convert_methrix(methrix_data)


gr <-  GenomicRanges::GRanges(paste0(rowData(methrix_data)[1,1], ":", rowData(methrix_data)[1,2], "-", rowData(methrix_data)[3,2]+1))
gr_df <- as.data.frame(gr)
gr_df2 <- gr_df
colnames(gr_df2)[1] <- "chr"


test_that("selected by regions", {
  expect_equal(region_filter(methrix_data, regions = gr), methrix_data[-(1:3),])
  expect_equal(region_filter(methrix_data, regions = gr_df2), methrix_data[-(1:3),])
})


test_that("Wrong input", {
  expect_error(region_filter("not methrix data", regions = gr))
  expect_error(region_filter(methrix_data, regions = "not granges"))
  expect_warning(region_filter(methrix_data, regions = gr_df))
})

test_that("selected by regions, HDF5", {
  expect_equal(region_filter(m2, regions = gr), m2[-(1:3),])
  expect_equal(region_filter(m2, regions = gr_df2), m2[-(1:3),])
})

