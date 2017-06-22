# library("devtools") ; library("roxygen2") ; library("fpca")
roxygenise() ; install(".")
library("fcomplete")
all.data = read.csv("/media/lukasz/DATA/alldata.csv")
all.data = all.data[(all.data$bmi < 100) & (all.data$height < 150) & (all.data$age) < 20,]
tb = table(all.data$Patient_ID)
#tb = tb[tb >= 3]
all.data = all.data[all.data$Patient_ID %in% names(tb),]

data = all.data[all.data$side=="L",c("Patient_ID","age","cadence")]
data = na.omit(data)
data = data[data$age < 20,]
dgrid = 100

mn = loess(data[,3] ~ data$age, intrapolate=FALSE)
grid = min(data$age) + (max(data$age) - min(data$age)) * (0:99)/100
pred = predict(mn, grid) #data.frame(age=grid))
#lines(grid,pred)

wide = fc.long2wide(data$Patient_ID, data$age, data[,3],bins = dgrid)
wide = t(t(wide) - pred)

long = fc.wide2long(wide)
plot(long$time,long$value)

# test:
# data[data$Patient_ID == 4030,2:3] - long[long$id == 4030,2:3]
# should give zeros in the second column and values smaller than (maxval - minval) / bins in the first column

smp = fc.sample(wide)
rm("fpca.model")

basis = fc.basis(5, "splines",norder = 4, dgrid = dgrid)

mean.impute = fc.mean(smp$train)
if (!("fpca.model" %in% ls())){
  long.train = fc.wide2long(smp$train)
  fpca.model = fc.fpca(long.train[,],d = 5, K = 2:4)
}
func.impute = functionalImpute(smp$train, fc.basis(5, "splines",norder = 4, dgrid = dgrid), maxIter = 1e4, thresh= 1e-5, lambda = fpca.model$sigma.est * c(0.1, 0.5, 0.7, 1, 1.5, 2) )

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

ensamble = (fpca.model$fit + func.impute$fit)/2
matplot(t(ensamble[smp$test.rows[20:25],]),t='l')
matplot(t(smp$train[smp$test.rows[20:25],]),t='p',pch = 'X',add=T)
matplot(t(smp$test[smp$test.rows[20:25],]),t='p',pch = 'O',add=T)

m1 = sqrt(mean((smp$test - func.impute$fit)[smp$test.mask]**2))
m0 = sqrt(mean((smp$test - mean.impute)[smp$test.mask]**2))
m2 = sqrt(mean((smp$test - fpca.model$fit)[smp$test.mask]**2))
m3 = sqrt(mean((smp$test - ensamble )[smp$test.mask]**2))

cat("mean impute:",m0,"\nfpca:",m2,"\nours:",m1,"\nensamble",m3)
