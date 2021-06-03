data('methrix_data')
m2 <- convert_methrix(methrix_data)

#' subset_methrix(methrix_data, contigs = 'chr21')
gr <-  GenomicRanges::GRanges(paste0(rowData(methrix_data)[1,1], ":", rowData(methrix_data)[1,2], "-", rowData(methrix_data)[3,2]+1))
gr_df <- as.data.frame(gr)
gr_df2 <- gr_df
colnames(gr_df2)[1] <- "chr"
samples <- rownames(colData(methrix_data))[c(1,3)]

test_that("selected by regions", {
  expect_equal(object = subset_methrix(methrix_data, regions = gr), expected = methrix_data[1:3,])
  expect_equal(subset_methrix(methrix_data, regions = gr_df2), methrix_data[1:3,])
  expect_equal(subset_methrix(methrix_data, regions = gr, samples = samples), methrix_data[1:3,c(1,3)])
  expect_equal(subset_methrix(methrix_data, contigs = "chr21"), methrix_data[which(rowData(methrix_data)$chr=="chr21"),])
})


test_that("Wrong input", {
  expect_error(subset_methrix("not methrix data", regions = gr))
  expect_error(subset_methrix(methrix_data, regions = "not granges"))
  expect_error(subset_methrix(methrix_data, contigs = "chr"))
  expect_error(subset_methrix(methrix_data, samples = "non existing"))
  expect_warning(subset_methrix(methrix_data, regions = gr_df))
  })

test_that("selected by regions", {
  expect_equal(subset_methrix(m2, regions = gr), m2[1:3,])
  expect_equal(subset_methrix(m2, regions = gr_df2), m2[1:3,])
  expect_equal(subset_methrix(m2, regions = gr, samples = samples), m2[1:3,c(1,3)])
  expect_equal(subset_methrix(m2, contigs = "chr21"), m2[which(rowData(m2)$chr=="chr21"),])
})


test_that("Wrong input", {
  expect_error(subset_methrix("not methrix data", regions = gr))
  expect_error(subset_methrix(m2, regions = "not granges"))
  expect_error(subset_methrix(m2, contigs = "chr"))
  expect_error(subset_methrix(m2, samples = "non existing"))
  expect_warning(subset_methrix(m2, regions = gr_df))
})
