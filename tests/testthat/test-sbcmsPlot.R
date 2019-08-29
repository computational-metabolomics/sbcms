context ("Test sbcmsPlot function")

test_that ("sbcmsPlot returns ggplot output", {

  classes <- sbcdata$class
  batch <- sbcdata$batch
  order <- c(1:nrow(sbcdata$data))
  data <- t(sbcdata$data[, 1:10])

  out <- QCRSC(df = data, order = order, batch = batch, classes = classes,
    spar = 0, minQC = 4)

  plots <- sbcmsPlot(df = data, corrected_df = out, 
    classes, batch, output=NULL, indexes=c(1,5,7))
  
  testthat::expect_s3_class(plots[[1]], "ggplot")
  # Empty plots are getting removed from the list object
  testthat::expect_true(length(plots) == 3)

})

test_that ("sbcmsPlot returns pdf output", {
  
  classes <- sbcdata$class
  batch <- sbcdata$batch
  order <- c(1:nrow(sbcdata$data))
  data <- t(sbcdata$data[, 1:2])
  
  out <- QCRSC(df = data, order = order, batch = batch, classes = classes,
    spar = 0, minQC = 4)
  pdf_output <- tempfile()
  testthat::expect_warning(sbcmsPlot (df = data, corrected_df = out, classes, 
    batch, output=pdf_output))
  
  testthat::expect_true(file.exists(pdf_output))
  testthat::expect_true(file.remove(pdf_output))
})

test_that ("sbcmsPlot returns output for the first 100 features if
  indexes=NULL", {
  
  classes <- sbcdata$class
  batch <- sbcdata$batch
  order <- c(1:nrow(sbcdata$data))
  data <- t(sbcdata$data[, 1:101])
  
  out <- QCRSC(df = data, order = order, batch = batch, classes = classes,
               spar = 0, minQC = 4)
  
  plots <- sbcmsPlot (df = data, corrected_df = out, classes, batch, output=NULL)
  
  testthat::expect_true(length(plots) == 100)
})
