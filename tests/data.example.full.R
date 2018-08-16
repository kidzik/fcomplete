library("devtools")
# library("roxygen2")
# roxygenise()
library("fpca")
install(".")
library("fcomplete")
library("ggplot2")
library("parallel")
source("tests/plot.helpers.R")

################
# PREPARE DATA #
################
if (!("all.data" %in% ls())){
  all.data = read.csv("/home/kidzik/Dropbox/DATA/CP/alldata.csv")
  gait.cycles = t(read.csv("/home/kidzik/Dropbox/DATA/CP/G_avg_CP.csv"))
  gdi = read.csv("/home/kidzik/Dropbox/DATA/CP/gdi.csv")
  pcas = prcomp(gait.cycles)
  all.data = cbind(all.data, pcas$x[,1:10])
  all.data = merge(all.data, gdi,by = c("Patient_ID","examid","side"))
#  all.data = read.csv("/home/lukasz/alldata.csv")
#  gait.cycles = t(read.csv("/home/lukasz/G_avg_CP.csv"))
}

all.data.subset = all.data[all.data$side == "L",]
# all.data.subset = all.data.subset[all.data.subset$age < 18,] # remove outliers
# all.data.subset = all.data.subset[all.data.subset$age > 5,] # remove outliers
all.data.subset = all.data.subset[all.data.subset$bmi > 5,] # remove outliers
all.data.subset = all.data.subset[all.data.subset$bmi < 25,] # remove outliers

cmeans = colMeans(all.data.subset[,300:309])
all.data.subset$score = sqrt(rowSums(t(t(all.data.subset[,300:309]) - cmeans)**2))

all.data.filtered = all.data.subset[,c("Patient_ID", "speed","cadence","age","bmi","height","KneeFlex_maxExtension","O2cost","GDI","PC1","PC2","PC3")]
all.data.filtered = all.data.filtered[order(all.data.filtered$Patient_ID, all.data.filtered$age),]
TBL = table(all.data.filtered$Patient_ID)
all.data.filtered = all.data.filtered[all.data.filtered$Patient_ID %in% names(TBL[TBL > 2]),]
all.data.filtered = na.omit(all.data.filtered)
all.data.filtered = all.data.filtered[!is.nan(all.data.filtered$GDI),]
pats = table(all.data.filtered$Patient_ID)
pats = names(pats[pats>=3])
all.data.filtered.sample = all.data.filtered[all.data.filtered$Patient_ID %in% pats,]

#################
# RUN CROSS-VAL #
#################
experiment.data = function(i)
{
  set.seed(i+20)
  # Sample data for testing

  var = "GDI"
  data = sample.long(all.data.filtered.sample, "Patient_ID", "age", var, ratio = 0.05, min.per.sbj = 3)

  # Set up parameters
  lambdas = list()
  lambdas[["GDI"]] = 5#seq(18,22,length.out = 5)
  lambdas[["bmi"]] = seq(1,2,length.out = 5)
  d = 6
  K = d

  # IMPUTE
  model.impute = fregression(as.formula(paste0(var,":age ~ 1 | Patient_ID")), data$train,
                             lambda= lambdas[[var]], thresh = 1e-10, maxIter = 10000,
                             method = "fimpute", final = "soft",
                             K=1, d=d, fold = 5)
  model.impute.fpcs = fregression(as.formula(paste0(var,":age ~ 1 | Patient_ID")), data$train, lambda= c(7.5), thresh = 1e-4, method = "fpcs", K=2:2, d=d)

  # lambdas[["bmi"]] = 0.25
  # X = list(train = data$train.matrix)
  # model.impute = fcomplete:::functionalMultiImpute.one(X, basis=model.impute.fpcs$basis, K=1, maxIter=10000, thresh=1e-10, lambda=lambdas[[var]])


#  model.mean = fregression(GDI:age ~ 1 | Patient_ID, data$train,  method = "mean")

  # REGRESSION
  # lambdas.reg = seq(0,1,length.out = 5)
  model.regression = fregression(GDI:age ~ GDI + O2cost + speed | Patient_ID, data$train,
                                method = "fimpute", thresh=1e-10, maxIter = 5000,
                                K.reg = 1,
                                lambda = lambdas[[var]],
                                lambda.reg = 1,
                                d=d,
                                K=1)

  errors = c(mean((model.regression$fit - data$test.matrix)**2, na.rm = TRUE),
    mean((model.impute$fit - data$test.matrix)**2, na.rm = TRUE),
    mean((model.impute.fpcs$fit - data$test.matrix)**2, na.rm = TRUE),
    mean((mean(data$train.matrix,na.rm=TRUE) - data$test.matrix)**2, na.rm = TRUE)
    )
  names(errors) = c("reg","impute","fPCA","mean")
  print(errors)
  list(errors = errors,
       model.fimp = model.impute,
       model.fpca = model.impute.fpcs,
       model.fslr = model.regression,
#       model.mean = model.mean,
       data = data)
}

models = mclapply(1:20, experiment.data, mc.cores = 4)
#models = lapply(1:1, data.experiment)

if (length(models))
  save(models, file="data-study.Rda")
if (!length(models))
  load("data-study.Rda")

exp.id = 1

# Basic diagnostics
ind = 1:2 + 50
obs = models[[exp.id]]$model.fimp$Y[ind,]
plot_preds(obs, models[[exp.id]]$model.fimp$fit[ind,], models[[exp.id]]$model.fslr$fit[ind,],
           filename="pred-data", title = "Mean", d=51)
p = 0.7
sm = models[[exp.id]]$model.fslr$fitI*p + (1-p)*models[[exp.id]]$model.fslr$fitR
sm = t(t(sm) + models[[exp.id]]$model.fslr$cmeans)

res = c()
for (i in 1:length(models)){
  if (is.null(attr(models[[i]],"class")))
    res = cbind(res, models[[i]]$errors)
}
rowMeans(res)
apply(res,1,sd)

# Compare predictions (for debugging)
exp.id = 2

models[[exp.id]]$model.fimp$fit
models[[exp.id]]$model.fpca$fit[1,]

mse.fimp = rowSums((models[[exp.id]]$data$test.matrix - models[[exp.id]]$model.fimp$fit)**2,na.rm = TRUE)
mse.fpca = rowSums((models[[exp.id]]$data$test.matrix - models[[exp.id]]$model.fpca$fit)**2,na.rm = TRUE)
test.points = mse.fpca > 1e-15
plot(mse.fimp[test.points])
plot(mse.fpca[test.points])
which(mse.fimp[test.points] > 30)

models[[exp.id]]$data$X[models[[exp.id]]$data$X$Patient_ID == 4822,]

# Summarize results
rownames(res) = c("SLR","SLI","fPCA","mean")
cbind(rowMeans(res),
apply(res,FUN=sd,1))
colnames(res) = paste("run",1:length(models))
ind = row.names(models[[i]]$data$test.matrix) %in% models[[i]]$data$X$Patient_ID[models[[i]]$data$test.ob[1:3]]

# Figure 7: Boxplot of results
library(tidyr)
tmp = t(res)
rownames(tmp) = 1:nrow(tmp)
tmp[,1:3] = 1 - t(t(tmp[,1:3]) / tmp[,4])
methodStats = gather(data.frame(tmp), method, varexp, SLR:fPCA, factor_key = FALSE)

pp = ggplot(methodStats, aes(x = method, y = varexp)) + paper.theme + labs(x="Method",y="MSE") +
  geom_boxplot()
pp
myggsave(filename=paste0("docs/plots/data-boxplot.pdf"), plot=pp, width = 10, height = 8)

# Lambdas
models[[exp.id]]$model.fimp$meta
models[[exp.id]]$model.fpca$selected_model

#############################
# OTHER PLOTS FOR THE PAPER #
#############################
library("ggplot2")
library("ggthemes")
library(RColorBrewer)

dd = all.data.filtered[,c("Patient_ID","age","bmi","GDI")]
dd$Patient_ID = as.factor(dd$Patient_ID)
dd = dd[dd$bmi > 10,]

# Figure 1: BMI over time
pp = ggplot(aes(x = age, y = bmi, color = Patient_ID), data = dd[1:200,]) + ylab("BMI") +
  geom_point(size = 3) + theme_set(theme_grey(base_size = 26)) + theme(legend.position="none", panel.background = element_rect(fill = "white",linetype = 1,colour = "grey50",size = 1,)) +
  stat_function(fun = approxfun(lowess(dd$age,dd$bmi)), size = 1.5, colour = "#000000")+ scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0))
pp
ggsave("docs/plots/points.pdf",width=7,height=5)
pp + geom_line(size=0.7)
ggsave("docs/plots/grouped.pdf",width=7,height=5)

# # The palette with grey:
# # The palette with black:
# cbPalette <- c("#000000", "#A69F00", "#56B4E9", "#009E73", "#008442") #, "#0072B2", "#D55E00", "#CC79A7")
#
# # Figure 3: Compare PCs with singular vectors (GDI over time)
# gg = as.numeric(colnames(model.impute.fpcs$Y))
# pp = ggplot(aes(x = age, y = GDI, color = Patient_ID), data = dd[2:14,]) +
#   scale_fill_manual(values = cbPalette) + scale_colour_manual(values = cbPalette) +
#   ylab("1st component") +
#   geom_point(size = 3) + theme_set(theme_grey(base_size = 26)) + theme(legend.position="none", panel.background = element_rect(fill = "white",linetype = 1,colour = "grey50",size = 1,)) +
#   stat_function(fun = approxfun(lowess(dd$age,dd$GDI)), size = 2.5, colour = "#000000")+ scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0)) + xlim(5,13)
# pp + geom_line(size = 0.7) +
#   stat_function(fun = approxfun(gg, model.impute.fpcs$fit[5,]) , color = cbPalette[5]) +
#   stat_function(fun = approxfun(gg, model.impute.fpcs$fit[2,]) , color = cbPalette[2]) +
#   stat_function(fun = approxfun(gg, model.impute.fpcs$fit[3,]) , color = cbPalette[3])  +
#   stat_function(fun = approxfun(gg, model.impute.fpcs$fit[4,]) , color = cbPalette[4])
