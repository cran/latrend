context('lcModel implementation')

setClass('lcModelTest', contains = 'lcModel')
testModel = latrend(lcMethodTestKML(), data = testLongData)
class(testModel) = 'lcModelTest'

predClusFun = function(object, newdata = NULL, cluster, ...) {
  rep(NaN, nrow(newdata))
}

predFun = function(object, newdata, ...) {
  pred = matrix(NaN, nrow = nrow(newdata), ncol = nClusters(object))
  transformPredict(pred = pred, model = object, newdata = newdata)
}

# including this test results in error for predict() and fitted() in later tests. No clue why.
# test_that('no predict funs', {
#   expect_error(predict(testModel, newdata=data.frame(Assessment=1)))
#   expect_error(predictForCluster(testModel, newdata=data.frame(Assessment=1), cluster = 'A'))
# })

setMethod('predictForCluster', signature('lcModelTest'), definition = predClusFun)

test_that('default predict.lcModel', {
  dfpred = predict(testModel, newdata=data.frame(Assessment=1))

  expect_is(dfpred, 'list')
  expect_is(dfpred$A$Fit, 'numeric')
  expect_equivalent(nrow(dfpred$A), 1)

  # removeMethod('predictForCluster', 'lcModelTest')
})

# NOTE: disabled until there is a way to unregister an S3 method
# test_that('default predictForCluster', {
#   .S3method('predict', 'lcModelTest', predFun)
#   pred = predictForCluster(testModel, newdata=data.frame(Assessment=c(1,2)), cluster = 'A')
#   expect_is(pred, 'numeric')
#   expect_equ
#   .S3method('predict', 'lcModelTest', predict.lcModel)
# })


test_that('default fitted', {
  # setMethod('predictForCluster', signature('lcModelTest'), predClusFun)

  suppressWarnings({
    expect_is(fitted(testModel), 'numeric')
  })

  # removeMethod('predictForCluster', 'lcModelTest')
})

test_that('name', {
  name = getName(testModel)
  expect_equal(name, getName(lcMethodTestKML()))
})

test_that('shortname', {
  name = getShortName(testModel)
  expect_equal(name, getShortName(lcMethodTestKML()))
})