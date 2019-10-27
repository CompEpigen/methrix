
library(GenomicRanges)

gr <- GenomicRanges::GRanges(paste0(methrix_data@elementMetadata[1,1], ":", methrix_data@elementMetadata[1,2], "-", methrix_data@elementMetadata[3,2]))
data("methrix_data")


M_mean <- apply(get_matrix(methrix_data[1:3,], type = "M", add_loci = FALSE), 2, mean, na.rm=TRUE)
C_mean <- apply(get_matrix(methrix_data[1:3,], type = "C", add_loci = FALSE), 2, mean, na.rm=TRUE)


M_min <- apply(get_matrix(methrix_data[1:3,], type = "M", add_loci = FALSE), 2, min, na.rm=TRUE)
C_min <- apply(get_matrix(methrix_data[1:3,], type = "C", add_loci = FALSE), 2, min, na.rm=TRUE)


M_max <- apply(get_matrix(methrix_data[1:3,], type = "M", add_loci = FALSE), 2, max, na.rm=TRUE)
C_max <- apply(get_matrix(methrix_data[1:3,], type = "C", add_loci = FALSE), 2, max, na.rm=TRUE)


test_that("selected regions", {
  expect_equal(as.numeric(get_region_summary(methrix_data,
                                                regions=gr, type = "M", how = "mean", overlap_type = "any")[1,-(1:4)]), as.numeric(M_mean))
  expect_equal(as.numeric(get_region_summary(methrix_data,
                                                regions=gr, type = "C", how = "mean", overlap_type = "any")[1,-(1:4)]), as.numeric(C_mean))

  expect_equal(as.numeric(get_region_summary(methrix_data,
                                                regions=gr, type = "M", how = "min", overlap_type = "any")[1,-(1:4)]), as.numeric(M_min))
  expect_equal(as.numeric(get_region_summary(methrix_data,
                                                regions=gr, type = "C", how = "min", overlap_type = "any")[1,-(1:4)]), as.numeric(C_min))

  expect_equal(as.numeric(get_region_summary(methrix_data,
                                                regions=gr, type = "M", how = "max", overlap_type = "any")[1,-(1:4)]), as.numeric(M_max))
  expect_equal(as.numeric(get_region_summary(methrix_data,
                                                regions=gr, type = "C", how = "max", overlap_type = "any")[1,-(1:4)]), as.numeric(C_max))
})



test_that("Using within", {
  expect_false(all(as.numeric(get_region_summary(methrix_data,
                                                regions=gr, type = "M", how = "mean", overlap_type = "within")[1,-(1:4)])==as.numeric(M_mean)))
  expect_false(all(as.numeric(get_region_summary(methrix_data,
                                                regions=gr, type = "C", how = "mean", overlap_type = "within")[1,-(1:4)])==as.numeric(C_mean)))
})

test_that("Wrong input", {
  expect_error(get_region_summary(methrix_data, regions=gr, type = "T", how = "mean", overlap_type = "within"))
  expect_error(get_region_summary(methrix_data, regions="not regions", type = "M", how = "mean", overlap_type = "within"))
  expect_error(get_region_summary(methrix_data, regions=gr, type = "M", how = "not existing", overlap_type = "within"))
  expect_error(get_region_summary(methrix_data, regions=gr, type = "M", how = "mean", overlap_type = "not existing"))
  expect_error(get_region_summary("not methrix data", regions=gr, type = "M", how = "mean", overlap_type = "within"))
  })


