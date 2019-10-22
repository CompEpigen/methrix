data('methrix_data')
m2 <- convert_methrix(methrix_data)

#' subset_methrix(methrix_data, contigs = 'chr21')


gr <-  GRanges(paste0(methrix_data@elementMetadata[1,1], ":", methrix_data@elementMetadata[1,2], "-", methrix_data@elementMetadata[3,2]+1))
gr_df <- as.data.frame(gr)
gr_df2 <- gr_df
colnames(gr_df2)[1] <- "chr"
samples <- rownames(methrix_data@colData)[c(1,3)]

test_that("selected by regions", {
  expect_identical(subset_methrix(methrix_data, regions = gr), methrix_data[1:3,])
  expect_identical(subset_methrix(methrix_data, regions = gr_df2), methrix_data[1:3,])
  expect_identical(subset_methrix(methrix_data, regions = gr, samples = samples), methrix_data[1:3,c(1,3)])
  expect_identical(subset_methrix(methrix_data, contigs = "chr21"), methrix_data[which(methrix_data@elementMetadata$chr=="chr21"),])
})


test_that("Wrong input", {
  expect_error(subset_methrix("not methrix data", regions = gr))
  expect_error(subset_methrix(methrix_data, regions = "not granges"))
  expect_error(subset_methrix(methrix_data, contigs = "chr"))
  expect_error(subset_methrix(methrix_data, samples = "non existing"))
  expect_warning(subset_methrix(methrix_data, regions = gr_df))
  })

test_that("selected by regions", {
  expect_identical(subset_methrix(m2, regions = gr), m2[1:3,])
  expect_identical(subset_methrix(m2, regions = gr_df2), m2[1:3,])
  expect_identical(subset_methrix(m2, regions = gr, samples = samples), m2[1:3,c(1,3)])
  expect_identical(subset_methrix(m2, contigs = "chr21"), m2[which(m2@elementMetadata$chr=="chr21"),])
})


test_that("Wrong input", {
  expect_error(subset_methrix("not methrix data", regions = gr))
  expect_error(subset_methrix(m2, regions = "not granges"))
  expect_error(subset_methrix(m2, contigs = "chr"))
  expect_error(subset_methrix(m2, samples = "non existing"))
  expect_warning(subset_methrix(m2, regions = gr_df))
})
