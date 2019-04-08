context ("test-splineSmoother function")

test_that ("splineSmoother returns expected output", {

  y_qc <- sbcdata$data[sbcdata$class=="QC" ,1]
  x_qc <- which(sbcdata$class=="QC")
  x_sample <- which(sbcdata$class!="QC")
  out <- splineSmoother(x=x_qc, y=y_qc, newX=x_sample, log=TRUE, spar=0, a=1)
  expect_equal (out, testData$splineSmoother)

})
