roxygenise() ; install(".")
library("fcomplete")
#all.data = read.csv("/media/lukasz/TOSHIBA EXT/alldata.csv")
data = all.data[all.data$side=="L",c("Patient_ID","age","bmi")]
data = na.omit(data)
data = data[data$age < 30,]
dgrid = 50

wide = fc.long2wide(data$Patient_ID, data$age, data$bmi,bins = dgrid)
#long = fc.wide2long(wide)

# test:
# data[data$Patient_ID == 4030,2:3] - long[long$id == 4030,2:3]
# should give zeros in the second column and values smaller than (maxval - minval) / bins in the first column

smp = fc.sample(wide)

basis = fc.basis(11,type = "splines",norder = 4, dgrid = dgrid)
R = functionalImpute(smp$train, basis, K=2)

# Plot
matplot(t(R$fit[smp$test.rows[20:25],]),t='l')
matplot(t(smp$test[smp$test.rows[20:25],]),t='p',pch = 'X',add=T)

m1 = sqrt(mean((smp$test - R$fit)[smp$test.mask]**2))

m0 = sqrt(mean((smp$test - fc.mean(smp$train))[smp$test.mask]**2))
cat(paste("Variance explained:",(1 - m1/m0)*100,"%"))

long.train = fc.wide2long(smp$train)
#long.test = fc.wide2long(smp$test)

model = fc.fpca(long.train)
