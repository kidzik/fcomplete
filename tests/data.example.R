library("devtools") ; library("roxygen2") ; library("fpca")
roxygenise() ; install(".")
library("fcomplete")
all.data = read.csv("/media/lukasz/DATA/alldata.csv")
all.data = all.data[(all.data$bmi < 100) & (all.data$height < 150) & (all.data$age) < 20,]
tb = table(all.data$Patient_ID)
#tb = tb[tb >= 3]
all.data = all.data[all.data$Patient_ID %in% names(tb),]

data = all.data[all.data$side=="L",c("Patient_ID","age","speed")]
data = na.omit(data)
data = data[(data$age < 15) & (data$age > 4),]
dgrid = 50

mn = loess(data[,3] ~ data$age, intrapolate=FALSE)
grid = min(data$age) + (max(data$age) - min(data$age)) * (0:(dgrid-1))/(dgrid-1)
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

nsplines = 6
#basis = fc.basis(nsplines, "splines",norder = 4, dgrid = dgrid)

mean.impute = fc.mean(smp$train)
if (!("fpca.model" %in% ls())){
  long.train = fc.wide2long(smp$train)
  fpca.model = fc.fpca(long.train[,],d = nsplines, K = 2:4, grid)
}
sigma.factors = c(0.1, 0.5, 0.7, 1, 1.5, 2, 2.5)
func.impute = functionalImpute(smp$train,
                               fc.basis(nsplines, "splines",norder = 4, dgrid = dgrid),
                               maxIter = 1e5, thresh= 1e-4,
                               lambda = fpca.model$sigma.est * fpca.model$sigma.est * sigma.factors)

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

m0 = sqrt(mean((smp$test - mean.impute)[smp$test.mask]**2))
m1 = sqrt(mean((smp$test - func.impute$fit)[smp$test.mask]**2))
m2 = sqrt(mean((smp$test - fpca.model$fit)[smp$test.mask]**2))
m3 = sqrt(mean((smp$test - ensamble )[smp$test.mask]**2))

cat("mean impute:\t",m0,"\nours:\t\t",m1,"\nfpca:\t\t",m2,"\nensamble:\t",m3)
