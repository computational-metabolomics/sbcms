context ("Test QCRSC function")

test_that ("QC-RSC returns expected output", {

  classes <- sbcdata$class
  batch <- sbcdata$batch
  order <- c(1:nrow(sbcdata$data))
  data <- t(sbcdata$data[, 1:20])

  out <- QCRSC(df = data, order = order, batch = batch, classes = classes,
                  spar = 0, minQC = 4)

  expect_equal (out, testData$QCRSC)

})
