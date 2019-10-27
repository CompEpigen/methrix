library(methrix)
library(GenomicRanges)
data("methrix_data")

m1 <- methrix_data
m2 <- convert_methrix(methrix_data)

dt <- as.data.frame(cbind(m1@elementMetadata, m1@assays$data$beta))
data.table::setDT(x = dt)
dt2 <- assays(m1)$beta

dt_c <- setDT(as.data.frame(cbind(m1@elementMetadata, m1@assays$data$cov)))
dt2_c <- assays(m1)$cov


test_that("Check output", {
  expect_equivalent(get_matrix(m = m1, type = "M",add_loci = FALSE, in_granges = FALSE), dt2)
  expect_equivalent(get_matrix(m = m1, type = "M",add_loci = TRUE, in_granges = FALSE), dt)
  expect_equivalent(get_matrix(m = m1, type = "C",add_loci = FALSE, in_granges = FALSE), dt2_c)
  expect_equivalent(get_matrix(m = m1, type = "C",add_loci = TRUE, in_granges = FALSE), dt_c)
})
dt_gr <- as.data.frame(cbind(m1@elementMetadata, m1@assays$data$beta))
dt_gr$end <- dt_gr$start+2
dt_gr <- GenomicRanges::makeGRangesFromDataFrame(dt_gr, keep.extra.columns = TRUE)

dt_c_gr <- as.data.frame(cbind(m1@elementMetadata, m1@assays$data$cov))
dt_c_gr$end <- dt_c_gr$start+2
dt_c_gr <- GenomicRanges::makeGRangesFromDataFrame(dt_c_gr, keep.extra.columns = TRUE)


test_that("Check GRanges output", {
  expect_equivalent(get_matrix(m = m1, type = "M",add_loci = TRUE, in_granges = TRUE), dt_gr)
  expect_equivalent(get_matrix(m = m1, type = "C",add_loci = TRUE, in_granges = TRUE), dt_c_gr)
})


test_that("Check HDF5 output", {
  expect_is(get_matrix(m = m2, type = "M",add_loci = FALSE, in_granges = FALSE), "DelayedMatrix")
  expect_is(get_matrix(m = m2, type = "C",add_loci = FALSE, in_granges = FALSE), "DelayedMatrix")
  })

test_that("Check wrong input", {
  expect_error(get_matrix(m = m1, type = "T", add_loci = FALSE, in_granges = FALSE))
  expect_error(get_matrix(m = "not methrix", type = "T", add_loci = FALSE, in_granges = FALSE))
  expect_error(get_matrix(m = "not methrix", type = "T", add_loci = FALSE, in_granges = TRUE))
  expect_warning(get_matrix(m = m1, type = "M", add_loci = FALSE, in_granges = TRUE))
})
