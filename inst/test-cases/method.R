# Testing the core functionality of the lcMethod class
####

mTest = make.lcMethod(id = 'Test.Id', time = 'Test.Time', response = 'Test.Response')

test('getName.nzchar', nchar(getName(mTest)) > 0)
test('getShortName.nzchar', nchar(getShortName(mTest)) > 0)
test('getLabel.ischar', is.character(getLabel(mTest)))
test('idVariable', idVariable(mTest), 'Test.Id')
test('timeVariable', timeVariable(mTest), 'Test.Time')
test('responseVariable', responseVariable(mTest), 'Test.Response')
test('getArgumentDefaults', is.list(getArgumentDefaults(mTest)))
test('getArgumentExclusions', is.character(getArgumentExclusions(mTest)))
test('getCall', is.call(getCall(mTest)))
test('getCall.fun', as.character(getCall(mTest)[[1]]) == class(mTest)[1])
test('getCitation', inherits(getCitation(mTest), 'citation'))

mTest3 = make.lcMethod(id = 'ID', time = 'TIME', response = 'Y', nClusters = 3L)
test('idVariable.3', idVariable(mTest3), 'ID')
test('timeVariable.3', timeVariable(mTest3), 'TIME')
test('responseVariable.3', responseVariable(mTest3), 'Y')
test('nClusters.3', mTest3$nClusters, 3L)
