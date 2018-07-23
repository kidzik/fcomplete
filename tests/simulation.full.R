# Regenerate and reinstall fimpute
library("roxygen2") ; roxygenize()
library("devtools")
devtools::install(".")
library("fcomplete")
library("ggplot2")
library("latex2exp")

res = list()
nexp = 1
dgrid = 31
d = 7

experiment.sim = function(exp.id){
  # SIMULATE DATA
  set.seed(exp.id)
  simulation = fsimulate(dgrid = dgrid,clear = 0.9, n = 100, noise.mag = 0.05, d = d, K = 1)
  data = simulation$data
  ftrue = simulation$ftrue
  K = simulation$params$K

  # TUNING PARAMS
  lambdas.pca = seq(0,2,length.out = 20)
  lambdas.reg = seq(0,2,length.out = 10)

  model.mean = fregression(Y:time ~ 1 | id, data, method = "mean", bins = dgrid)
  # model.fpca = fregression(Y:time ~ 1 | id, data, lambda = 0, K = 2:d, thresh = 1e-7, method = "fpcs", bins = dgrid)
  # model.fimp = fregression(Y:time ~ 1 | id, data, lambda = lambdas.pca, K = d, thresh = 0, final = "soft", maxIter = 1000, fold = 5, cv.ratio = 0.05, bins = dgrid)
  # model.fcmp = fregression(0:time ~ Y + X1 + X2 | id, data, lambda = lambdas.pca, final = "soft", bins = dgrid)
  # model.freg = fregression(Y:time ~ X1 + X2 | id, data, lambda = lambdas.reg, thresh = 1e-4, lambda.reg = 0.1 * 1:20, method = "fimpute", bins = dgrid)
  model.fslr = fregression(Y:time ~ Y + X1 + X2 | id, data, lambda = lambdas.reg, thresh = 1e-4, lambda.reg = 0.1 * 1:20, method = "fimpute", bins = dgrid)

  # REPORT RESULTS
  errors = c(
    mean((ftrue - mean(data$Y))**2),
    # mean((ftrue - model.fpca$fit)**2),
    # mean((ftrue - model.fimp$fit)**2),
    # mean((ftrue - model.fcmp$fit)**2),
    mean((ftrue - model.fslr$fit)**2)
  )
  tbl.true = cbind(
    errors,
    100*(1-errors/errors[1])
  )
  colnames(tbl.true) = c("MSE","% expl")
  rownames(tbl.true) = c("mean",
#                         "fpca","fimpute","fcompress",
                         "regression")
  print(tbl.true)

  # SAVE EXPERIMENT RESULTS
  res = list()
  res$tbl = tbl.true
  res$errors = errors
  res$mean = model.mean
  # res$fimp = model.fimp
  # res$fpca = model.fpca
  res$fslr = model.fslr
  res$simulation = simulation
  res
}

res = lapply(1:5, experiment.sim)

#save(res,file = "sim-study.Rda")
#load("sim-study.Rda")
rowMeans(boxplot.data[,-c(4,6)])

exp.id = 1
#1 - res[[exp.id]]$errors /res[[exp.id]]$errors[1]

## Post-process results
boxplot.data = c()
for (i in 1:length(res)){
  boxplot.data = cbind(boxplot.data, res[[i]]$errors)
}
rownames(boxplot.data) = rownames(res[[exp.id]]$tbl)
joint.tbl = rowMeans(boxplot.data)

# Sanity check for the regression result
# p = 0.2
# ps = 0:200 / 100
# er = sapply(ps, function(p){
# fit = (model.fslr$fitI * (2-p) + model.fslr$fitR * p)
# mean((simulation$ftrue - fit)**2,na.rm = TRUE)
# })
# plot(ps,er)
#
# fit = (model.fslr$fitI + model.fslr$fitR)
# mean((simulation$fobs - model.fslr$fit)**2,na.rm = TRUE)
# mean((ftrue - fit)**2)

################
# PLOT RESULTS #
################
source("tests/plot.helpers.R")

# Figure 3:
exp.show = 8
res[[exp.show]]$tbl
#matplot(t(-res[[exp.show]]$fpca$v)[,1:min(3,nrow(res[[exp.show]]$fpca$v))],t='l',lwd = 4)

scalar = sqrt(rowSums(res[[exp.show]]$fpca$v**2)[1])
vlist = list(res[[exp.show]]$fpca$v / scalar,
             res[[exp.show]]$fimp$v)
vlist[[2]][2:3,] = vlist[[2]][2:3,]*(-1)
ttl = c("Components of fPCA","Components of SFI")
names = c("fpca","fimp")

for (j in 1:length(vlist)){
  cols = gg_color_hue(3)
  v = vlist[[j]]
  pp = ggplot(data.frame(x=c(0, 1)), aes(x)) + paper.theme +
    xlim(0,1) + ylim(-0.4,0.5) + labs(x = "time", y = "value", title=ttl[j])
  for (i in 1:min(3,nrow(v))){
    df = data.frame(x=0:30/30,y=v[i,])
    pp = pp + geom_segment(data=df, aes(x=x,y=y,xend=dplyr::lead(x),yend=dplyr::lead(y)),
                           color=cols[i], size=1.5, linetype=i)
  }
  print(pp)
  myggsave(filename=paste0("docs/plots/components-",names[[j]],".pdf"), plot=pp)
}

# Figure 4:
ind = 1:3 + 10
# plot_preds(res[[exp.show]]$simulation$fobs[ind,], res[[exp.show]]$simulation$ftrue[ind,], res[[exp.show]]$mean$fit[ind,],
#            filename="pred-mean", title = "Mean")
ylim = c(-0.5,1.5)
plot_preds(res[[exp.show]]$simulation$fobs[ind,], res[[exp.show]]$simulation$ftrue[ind,], NULL,
           filename="pred-mean", title = "True curves & observations",ylim=ylim)
plot_preds(res[[exp.show]]$simulation$fobs[ind,], res[[exp.show]]$simulation$ftrue[ind,], res[[exp.show]]$fpca$fit[ind,],
           filename="pred-fpca", title = "Functional PCA",ylim=ylim)
plot_preds(res[[exp.show]]$simulation$fobs[ind,], res[[exp.show]]$simulation$ftrue[ind,], res[[exp.show]]$fimp$fit[ind,],
           filename="pred-fimp", title = "Sparse Functional Impute",ylim=ylim)
plot_preds(res[[exp.show]]$simulation$fobs[ind,], res[[exp.show]]$simulation$ftrue[ind,], res[[exp.show]]$fslr$fit[ind,],
           filename="pred-freg", title = "Sparse Functional Regression",ylim=ylim)

# Figure 5: Sample error by cross-val
res[[exp.show]]$fimp$meta$cv.K = round(res[[exp.show]]$fimp$meta$cv.K)
pp = ggplot(data = res[[exp.show]]$fimp$meta, aes(x=lambda, y=cv.err)) + paper.theme +
  geom_line(size=1) +
  geom_point(size=4) +
  ggtitle("Error by cross-validation") +
  labs(x = TeX("$\\lambda$"), y = "error")
print(pp)
myggsave(filename=paste0("docs/plots/error-of-lambda.pdf"), plot=pp)

pp = ggplot(data = res[[exp.show]]$fimp$meta, aes(x=lambda, y=cv.K)) + paper.theme +
  geom_point(size=4) +
  ggtitle(TeX("Dimension K corresponding to $\\lambda$")) +
  labs(x = TeX("$\\lambda$"), y = "K")
print(pp)
myggsave(filename=paste0("docs/plots/K-of-lambda.pdf"), plot=pp)

# Figure 6: Example curves from the simulation
cols = gg_color_hue(3)

for (showLines in 0:1){
  pp = ggplot(data.frame(x=c(0, 1)), aes(x)) + paper.theme + ylim(-6,6) +
    xlim(0,1) + labs(x = "time", y = "value")
  for (i in 1:3){
    df = data.frame(x=0:30/30,y=simulation$ftrue[i,],yobs=res[[exp.show]]$simulation$fobs[i,])
    if (showLines)
      pp = pp + geom_segment(data=df, aes(x=x,y=y,xend=dplyr::lead(x),yend=dplyr::lead(y)),
                            color=cols[i], size=0.5)
    df = na.omit(df)
    if (showLines)
      pp = pp + geom_segment(data=df, aes(x=x,y=y,xend=x,yend=yobs),
                            color=cols[i], size=0.75, linetype="dotted")
    pp = pp + geom_point(data=df, aes(x=x,y=yobs),
                           color=cols[i], size=3)
  }
  print(pp)
}

# Figure 7: Boxplots for all methods
library(tidyr)
tmp = t(boxplot.data)[-c(4,6),-c(1,4)]
colMeans(tmp)
apply(tmp,2,sd)

rownames(tmp) = 1:nrow(tmp)
methodStats = gather(data.frame(tmp), method, varexp, fpca:regression, factor_key = FALSE)

pp = ggplot(methodStats, aes(x = method, y = varexp)) + paper.theme + labs(x="Method",y="Variance explained (%)") +
  geom_boxplot()
myggsave(filename=paste0("docs/plots/simulation-boxplot.pdf"), plot=pp, width = 12, height = 8)
pp
