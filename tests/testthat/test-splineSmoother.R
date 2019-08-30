context ("test-splineSmoother function")

test_that ("splineSmoother returns expected output", {

  y_qc <- sbcdata$data[sbcdata$class == "QC" & sbcdata$batch == 1 , 3]
  x_qc <- which(sbcdata$class=="QC" & sbcdata$batch == 1)
  x_sample <- which(sbcdata$class!="QC" & sbcdata$batch == 1)
  out <- sbcms:::splineSmoother(x=x_qc, y=y_qc, newX=x_sample,
    log=TRUE, spar=0, a=1)
  expect_equal (out, testData$splineSmoother)

})

test_that ("splineSmoother returns expected output and verbose message", {
  
  y_qc <- sbcdata$data[sbcdata$class == "QC" & sbcdata$batch == 1 , 3]
  x_qc <- which(sbcdata$class=="QC" & sbcdata$batch == 1)
  x_sample <- which(sbcdata$class!="QC" & sbcdata$batch == 1)
  out <- sbcms:::splineSmoother(x=x_qc, y=y_qc, newX=x_sample,
                                log=TRUE, spar=0, a=1, verbose=TRUE)
  expect_equal (out, testData$splineSmoother)
  
})