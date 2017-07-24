library("devtools") ; library("roxygen2") ; library("fpca")
roxygenise() ; install(".")
library("fcomplete")
all.data = read.csv("/media/lukasz/DATA/alldata.csv")
all.data = all.data[(all.data$bmi < 100) & (all.data$height < 150) & (all.data$age) < 20,]
tb = table(all.data$Patient_ID)
#tb = tb[tb >= 3]
all.data = all.data[all.data$Patient_ID %in% names(tb),]

data = all.data[all.data$side=="L",c("Patient_ID","age","bmi","speed")]
data = na.omit(data)
data = data[data$age < 20,]
dgrid = 100

wideX = fc.long2wide(data$Patient_ID, data$age, data[,3],bins = dgrid)
smp = fc.sample(wide)
func.impute.X = functionalImpute(smp$train, fc.basis(5, "splines",norder = 4, dgrid = dgrid), maxIter = 1e4, thresh= 1e-4, lambda = 1 )

wideY = fc.long2wide(data$Patient_ID, data$age, data[,4],bins = dgrid)
smp = fc.sample(wide)
func.impute.Y = functionalImpute(smp$train, fc.basis(5, "splines",norder = 4, dgrid = dgrid), maxIter = 1e4, thresh= 1e-4, lambda = 1 )

M = lm(func.impute.Y$U ~ func.impute.X$U)
summary(M)
