library("devtools")
# library("roxygen2")
# roxygenise()
library("fpca")
install(".")
library("fcomplete")
library("ggplot2")
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
all.data.subset = all.data.subset[all.data.subset$age < 15,] # remove outliers
all.data.subset = all.data.subset[all.data.subset$bmi > 5,] # remove outliers

cmeans = colMeans(all.data.subset[,300:309])
all.data.subset$score = sqrt(rowSums(t(t(all.data.subset[,300:309]) - cmeans)**2))

all.data.filtered = all.data.subset[,c("Patient_ID", "speed","cadence","age","bmi","height","KneeFlex_maxExtension","O2cost","GDI","PC1","PC2","PC3")]
all.data.filtered = all.data.filtered[order(all.data.filtered$Patient_ID, all.data.filtered$age),]
TBL = table(all.data.filtered$Patient_ID)
all.data.filtered = all.data.filtered[all.data.filtered$Patient_ID %in% names(TBL[TBL > 2]),]
all.data.filtered = na.omit(all.data.filtered)
all.data.filtered = all.data.filtered[!is.nan(all.data.filtered$GDI),]

#################
# RUN CROSS-VAL #
#################
data.experiment = function(i){
  set.seed(i)
  # Sample data for testing
  data = sample.long(all.data.filtered, "Patient_ID", "age", "GDI", ratio = 0.05)

  # Set up parameters
  lambdas = seq(10,50,length.out = 5)
  K = 6

  # IMPUTE
#  model.impute = fregression(GDI:age ~ 1 | Patient_ID, data$train, lambda= lambdas, thresh = 0, maxIter = 5000, method = "fimpute", K=K, d=K, fold = 3)
  model.impute.fpcs = fregression(GDI:age ~ 1 | Patient_ID, data$train, lambda= c(7.5), thresh = 1e-4, method = "fpcs", K=2:K, d=K)
#  model.mean = fregression(GDI:age ~ 1 | Patient_ID, data$train,  method = "mean")

  # REGRESSION
  lambdas.reg = seq(0.5,5,length.out = 10)
  #model.regression = fregression(GDI:age ~ GDI + bmi + O2cost | Patient_ID, data$train, method = "fimpute", K = K, thresh=1e-5, K.reg = K, lambda = lambdas, lambda.reg = lambdas, d=K)
  model.regression = fregression(GDI:age ~ GDI + bmi + O2cost | Patient_ID, data$train,
                                 method = "fimpute", K = 2, thresh=1e-5, K.reg = 2,
                                 lambda = lambdas.reg,
                                 lambda.reg = lambdas.reg/10, d=K)

  errors = c(mean((model.regression$fit - data$test.matrix)**2, na.rm = TRUE),
    #mean((model.impute$fit - data$test.matrix)**2, na.rm = TRUE),
    mean((model.impute.fpcs$fit - data$test.matrix)**2, na.rm = TRUE)
    #mean((model.mean$fit - data$test.matrix)**2, na.rm = TRUE)
    )
  names(errors) = c("regression")#,"impute","fPCA","mean")

  list(errors = errors,
#       model.fimp = model.impute,
       model.fpca = model.impute.fpcs,
       model.fslr = model.regression,
#       model.mean = model.mean,
       data = data)
}

library(parallel)
models = lapply(1:1, data.experiment)

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
  res = cbind(res, models[[i]]$errors)
}

# Summarize results
rownames(res) = c("regression","impute","fPCA","mean")
cbind(rowMeans(res),
apply(res,FUN=sd,1))
colnames(res) = paste("run",1:length(models))
ind = row.names(models[[i]]$data$test.matrix) %in% models[[i]]$data$X$Patient_ID[models[[i]]$data$test.ob[1:3]]

# Figure 7: Boxplot of results
library(tidyr)
tmp = t(res)
rownames(tmp) = 1:nrow(tmp)
methodStats = gather(data.frame(tmp), method, varexp, regression:fPCA, factor_key = FALSE)

pp = ggplot(methodStats, aes(x = method, y = varexp)) + paper.theme + labs(x="Method",y="MSE") +
  geom_boxplot()
pp
myggsave(filename=paste0("docs/plots/data-boxplot.pdf"), plot=pp, width = 5, height = 4)

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

# Figure 1: BMI over time
pp = ggplot(aes(x = age, y = bmi, color = Patient_ID), data = dd[1:200,]) + ylab("1st component") +
  geom_point(size = 3) + theme_set(theme_grey(base_size = 26)) + theme(legend.position="none", panel.background = element_rect(fill = "white",linetype = 1,colour = "grey50",size = 1,)) +
  stat_function(fun = approxfun(lowess(dd$age,dd$bmi)), size = 1.5, colour = "#000000")+ scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0))
pp
ggsave("docs/plots/points.pdf",width=10,height=7)
pp + geom_line(size=0.7)
ggsave("docs/plots/grouped.pdf",width=10,height=7)

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
