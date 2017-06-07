# library("devtools") ; library("roxygen2")
# library("fpca")
roxygenise() ; install(".")
library("fcomplete")
all.data = read.csv("/media/lukasz/TOSHIBA EXT/alldata.csv")
all.data = all.data[(all.data$bmi < 100) & (all.data$height < 150),]

data = all.data[all.data$side=="L",c("Patient_ID","age","speed")]
data = na.omit(data)
data = data[data$age < 20,]
dgrid = 100

mn = loess(data$speed ~ data$age, intrapolate=FALSE)
grid = min(data$age) + (max(data$age) - min(data$age)) * (0:99)/100
pred = predict(mn, grid) #data.frame(age=grid))
#lines(grid,pred)

wide = fc.long2wide(data$Patient_ID, data$age, data$speed,bins = dgrid)
wide = t(t(wide) - pred)

long = fc.wide2long(wide)
plot(long$time,long$value)

# test:
# data[data$Patient_ID == 4030,2:3] - long[long$id == 4030,2:3]
# should give zeros in the second column and values smaller than (maxval - minval) / bins in the first column

smp = fc.sample(wide)

basis = fc.basis(5, type = "splines",norder = 4, dgrid = dgrid)

roxygenise() ; install(".")
func.impute = functionalImpute(smp$train, basis, maxIter = 100, thresh= 0.0001, lambda = .001, K=1)
mean.impute = fc.mean(smp$train)
#long.train = fc.wide2long(smp$train)
#fpca.model = fc.fpca(long.train[,])

# Plot
matplot(t(mean.impute[smp$test.rows[20:25],]),t='l')
matplot(t(smp$train[smp$test.rows[20:25],]),t='p',pch = 'X',add=T)
matplot(t(smp$test[smp$test.rows[20:25],]),t='p',pch = 'O',add=T)

matplot(t(func.impute$fit[smp$test.rows[20:25],]),t='l')
matplot(t(smp$train[smp$test.rows[20:25],]),t='p',pch = 'X',add=T)
matplot(t(smp$test[smp$test.rows[20:25],]),t='p',pch = 'O',add=T)

matplot(t(fpca.model$fit[smp$test.rows[20:25],]),t='l')
matplot(t(smp$train[smp$test.rows[20:25],]),t='p',pch = 'X',add=T)
matplot(t(smp$test[smp$test.rows[20:25],]),t='p',pch = 'O',add=T)

m1 = sqrt(mean((smp$test - func.impute$fit)[smp$test.mask]**2))
m0 = sqrt(mean((smp$test - mean.impute)[smp$test.mask]**2))
m2 = sqrt(mean((smp$test - fpca.model$fit)[smp$test.mask]**2))

cat("mean impute:",m0,"\nfpca:",m2,"\nours:",m1)
