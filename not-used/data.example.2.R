library("devtools") ; library("roxygen2") ; library("fpca")
roxygenise() ; install(".")
library("fcomplete")
#all.data = read.csv("/media/lukasz/DATA/alldata.csv")
all.data = all.data[(all.data$bmi < 100) & (all.data$height < 150) & (all.data$age) < 20 & (all.data$leglen) < 80 & (all.data$leglen) > 0,]
tb = table(all.data$Patient_ID)
#tb = tb[tb >= 3]
all.data = all.data[all.data$Patient_ID %in% names(tb),]

data = all.data[all.data$side=="L",c("Patient_ID","age","PelvicTilt_mean","bmi","leglen")]
data[,3] = data[,3] / 17
data[,4] = data[,4] / 17
data[,5] = data[,5] / 60
data = na.omit(data)
data = data[(data$age < 15) & (data$age > 4),]
dgrid = 50

mn = loess(data[,3] ~ data$age, intrapolate=FALSE)
grid = min(data$age) + (max(data$age) - min(data$age)) * (0:(dgrid-1))/(dgrid-1)
#lines(grid,pred)

wide.X = fc.long2wide(data$Patient_ID, data$age, data[,3],bins = dgrid)
mn = loess(data[,3] ~ data$age, intrapolate=FALSE)
pred = predict(mn, grid) #data.frame(age=grid))
wide.X = t(t(wide.X) - pred)
wide.Y = fc.long2wide(data$Patient_ID, data$age, data[,4],bins = dgrid)
wide.Y = t(t(wide.Y))
wide.Z = fc.long2wide(data$Patient_ID, data$age, data[,5],bins = dgrid)
wide.Z = t(t(wide.Z))

smp.Y = fc.sample(wide.Y)

apply.mask = function(wide, mask){
  smp.X = list(test.mask = mask$test.mask)
  smp.X$train = wide
  smp.X$test = wide
  smp.X$test.rows = mask$test.rows
  nas = is.na(smp.X$train)
  smp.X$train[smp.X$test.mask] = NA
  smp.X$test[!smp.X$test.mask] = NA
  smp.X
}

smp.X = apply.mask(wide.X, smp.Y)
smp.Z = apply.mask(wide.Z, smp.Y)

nsplines = 6

rm("fpca.model")

mean.impute = fc.mean(smp.X$train)
if (!("fpca.model" %in% ls())){
  long.train = fc.wide2long(smp.X$train)
  fpca.model = fc.fpca(long.train[,],d = nsplines, K = 2, grid)
}

sigma.factors = c(0.1, 0.5, 0.7, 1, 1.5, 2, 2.5)

func.multi.impute = functionalMultiImpute.one(smp.X$train,
                                              basis = fc.basis(nsplines, "splines",norder = 4, dgrid = dgrid),
                                              maxIter = 1e2, thresh= 1e-4,
                                              lambda = 0,
                                              K = 2)
func.multi.impute.one = functionalImpute.one(smp.X$train,
                                              basis = fc.basis(nsplines, "splines",norder = 4, dgrid = dgrid),
                                              maxIter = 1e2, thresh= 1e-4,
                                              lambda = 0,
                                              K = 2)

func.impute.X = functionalImpute(smp.X$train,
                                     fc.basis(nsplines, "splines",norder = 4, dgrid = dgrid),
                                     maxIter = 1e4, thresh= 1e-3,
                                     lambda = fpca.model$sigma.est * fpca.model$sigma.est * sigma.factors)

func.impute.Y = functionalImpute(smp.Y$train,
                                          fc.basis(nsplines, "splines",norder = 4, dgrid = dgrid),
                                          maxIter = 1e2, thresh= 1e-4,
                                          lambda = 0,
                                          K = 2)




# Plot
matplot(t(mean.impute[smp.X$test.rows[20:25],]),t='l')
matplot(t(smp.X$train[smp.X$test.rows[20:25],]),t='p',pch = 'X',add=T)
matplot(t(smp.X$test[smp.X$test.rows[20:25],]),t='p',pch = 'O',add=T)

matplot(t(func.multi.impute[[1]][smp.X$test.rows[20:25],]),t='l')
matplot(t(smp.X$train[smp.X$test.rows[20:25],]),t='p',pch = 'X',add=T)
matplot(t(smp.X$test[smp.X$test.rows[20:25],]),t='p',pch = 'O',add=T)

matplot(t(func.multi.impute[[2]][smp.Y$test.rows[20:25],]),t='l')
matplot(t(smp.Y$train[smp.Y$test.rows[20:25],]),t='p',pch = 'X',add=T)
matplot(t(smp.Y$test[smp.Y$test.rows[20:25],]),t='p',pch = 'O',add=T)

matplot(t(func.impute.X$fit[smp.X$test.rows[20:25],]),t='l')
matplot(t(smp.X$train[smp.X$test.rows[20:25],]),t='p',pch = 'X',add=T)
matplot(t(smp.X$test[smp.X$test.rows[20:25],]),t='p',pch = 'O',add=T)

matplot(t(func.impute.Y$fit[smp.Y$test.rows[20:25],]),t='l')
matplot(t(smp.Y$train[smp.Y$test.rows[20:25],]),t='p',pch = 'X',add=T)
matplot(t(smp.Y$test[smp.Y$test.rows[20:25],]),t='p',pch = 'O',add=T)

#
matplot(t(fpca.model$fit[smp.X$test.rows[20:25],]),t='l')
matplot(t(smp.X$train[smp.X$test.rows[20:25],]),t='p',pch = 'X',add=T)
matplot(t(smp.X$test[smp.X$test.rows[20:25],]),t='p',pch = 'O',add=T)
#
# ensamble = (fpca.model$fit + func.impute$fit)/2
# matplot(t(ensamble[smp$test.rows[20:25],]),t='l')
# matplot(t(smp$train[smp$test.rows[20:25],]),t='p',pch = 'X',add=T)
# matplot(t(smp$test[smp$test.rows[20:25],]),t='p',pch = 'O',add=T)
#
m0 = sqrt(mean((smp.X$test - mean.impute)[smp.X$test.mask]**2))
m2 = sqrt(mean((smp.X$test - fpca.model$fit)[smp.X$test.mask]**2))

m.mX = sqrt(mean((smp.X$test - func.multi.impute[[2]])[smp.X$test.mask]**2))
m.mY = sqrt(mean((smp.Y$test - func.multi.impute[[1]])[smp.Y$test.mask]**2))

m.X = sqrt(mean((smp.X$test - func.impute.X$fit)[smp.X$test.mask]**2))
m.Y = sqrt(mean((smp.Y$test - func.impute.Y$fit)[smp.Y$test.mask]**2))

print(m0)
print(m2)
print(m.mX)
print(m.mY)
print(m.X)
print(m.Y)

# m3 = sqrt(mean((smp$test - ensamble )[smp$test.mask]**2))
#
# cat("mean impute:\t",m0,"\nours:\t\t",m1,"\nfpca:\t\t",m2,"\nensamble:\t",m3)
