context ("sbcWrapper function")

test_that ("sbcWrapper returns expected output", {

  classes <- sbcdata$class
  batch <- sbcdata$batch
  order <- c(1:nrow(sbcdata$data))

  Data <- t(sbcdata$data)
  qcData <- Data[,classes == "QC"]
  qc_batch <- batch[classes == "QC"]
  qc_order <- order[classes == "QC"]

  out <- sbcWrapper(id=4, qcData=qcData, order=order, qcOrder=qc_order, qcBatch=qc_batch,
                       log=TRUE, spar=0, batch=batch, minQC=4)

  expect_equal (out, testData$sbcWrapper)

})

test_that ("sbcWrapper returns expected output, with predefined spar value", {
  
  classes <- sbcdata$class
  batch <- sbcdata$batch
  order <- c(1:nrow(sbcdata$data))
  
  Data <- t(sbcdata$data)
  qcData <- Data[,classes == "QC"]
  qc_batch <- batch[classes == "QC"]
  qc_order <- order[classes == "QC"]
  
  out <- sbcWrapper(id=4, qcData=qcData, order=order, qcOrder=qc_order, qcBatch=qc_batch,
                    log=TRUE, spar=1, batch=batch, minQC=4)
  
  expect_equal (out, testData$sbcWrapper_spar_1)
  
})

test_that ("sbcWrapper returns NA if not enough data point are measured", {
  
  classes <- sbcdata$class
  batch <- sbcdata$batch
  order <- c(1:nrow(sbcdata$data))
  
  Data <- t(sbcdata$data)
  qcData <- Data[,classes == "QC"]
  qc_batch <- batch[classes == "QC"]
  qc_order <- order[classes == "QC"]
  
  qcData[4, -c(sample(38, 3))] <- NA
  
  out <- sbcWrapper(id=4, qcData=qcData, order=order, qcOrder=qc_order, qcBatch=qc_batch,
                    log=TRUE, spar=1, batch=batch, minQC=4)
  
  expect_true (all(is.na(out)))
  
})

test_that ("sbcWrapper returns error", {
  classes <- sbcdata$class
  batch <- sbcdata$batch
  order <- c(1:nrow(sbcdata$data))
  
  Data <- t(sbcdata$data)
  qcData <- Data[,classes == "QC"]
  qc_batch <- batch[classes == "QC"]
  qc_order <- order[classes == "QC"]
  
  dim(qcData) <- NULL
  
  expect_error(sbcWrapper(id=4, qcData=qcData, order=order, qcOrder=qc_order, qcBatch=qc_batch,
                    log=TRUE, spar=1, batch=batch, minQC=4))
})