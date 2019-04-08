context ("sbcWrapper function")

test_that ("sbcWrapper returns expected output", {

  classes <- sbcdata$class
  batch <- sbcdata$batch
  order <- c(1:nrow(sbcdata$data))

  Data <- t(sbcdata$data)
  qcData <- Data[,classes=="QC"]
  qc_batch <- batch[classes=="QC"]
  qc_order <- order[classes=="QC"]

  out <- sbcWrapper(id=4, qcData=qcData, order=order, qcOrder=qc_order, qcBatch=qc_batch,
                       log=TRUE, spar=0, batch=batch, minQC=4)

  expect_equal (out, testData$sbcWrapper)

})
