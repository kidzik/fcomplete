library("devtools") ; library("roxygen2") ; library("fpca")
roxygenise() ; install(".")
library("fcomplete")
all.data = read.csv("/media/lukasz/DATA/alldata.csv")
all.data = all.data[(all.data$bmi < 100) & (all.data$height < 150) & (all.data$age) < 20,] # remove outliers

# REGRESSION
model.regression = fregression(speed:age ~ leglen + bmi | Patient_ID, all.data)

# pid = 600
# p = all.data[all.data$Patient_ID == model$id[pid] & !is.na(all.data$Patient_ID),]
# row = model$fit[model$id == model$id[pid],]
# plot(model$grid, row,t='l', xlab = "age", ylab = "speed")
# points(p$age,p$speed)

# IMPUTE
model.impute = fregression(speed:age ~ 1 | Patient_ID, all.data)

# DIM REDUCTION
model.unsupervised = fregression(1:age ~ speed + leglen + bmi | Patient_ID, all.data)
