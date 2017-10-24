# Regenerate and reinstall fimpute
library("roxygen2") ; roxygenize()
library("devtools") ; devtools::install(".")
library("fcomplete")

#rm(list = ls())
res = list()
nexp = 1
dgrid = 51

for(exp.id in 1:nexp){

# SIMULATE DATA
set.seed(323 + exp.id)
simulation = fsimulate(dgrid = dgrid,clear = 0.95, n = 100,noise.mag = 0.03,)
data = simulation$data
ftrue = simulation$ftrue
K = 6 #simulation$params$K

# REGRESSION
model.mean = fregression(Y:time ~ 1 | id, data, method = "mean", bins = dgrid)
model.fpca = fregression(Y:time ~ 1 | id, data, lambda = 0, K = 2:K, thresh = 1e-7, method = "fpcs", bins = dgrid)

lambdas = seq(0,10,length.out = 20)
model.fimp = fregression(Y:time ~ 1 | id, data, lambda = lambdas, thresh = 0, final = "soft", maxIter = 100, fold = 5, cv.ratio = 0.05, K = K, bins = dgrid)
model.fcmp = fregression(0:time ~ Y + X1 + X2 | id, data, lambda = lambdas, K = K, final = "soft")
lambdas = seq(0,2,length.out = 10)
model.freg = fregression(Y:time ~ X1 + X2 | id, data, lambda = lambdas, thresh = 1e-5, lambda.reg = 0.1 * 5:20, method = "fimpute", K = K)

# REPORT RESULTS
idx = unique(data$id)
errors = c(
  mean((ftrue[idx,] - model.mean$fit)**2),
  mean((ftrue[idx,] - model.fpca$fit)**2),
  mean((ftrue[idx,] - model.fimp$fit)**2),
  mean((ftrue[idx,] - model.fcmp$fit)**2),
  mean((ftrue[idx,] - model.freg$fit)**2)
)
errors
model.fimp$meta

tbl.true = cbind(
  errors,
  100*(1-errors/errors[1])
)
colnames(tbl.true) = c("MSE","% expl")
rownames(tbl.true) = c("mean","fpca","fimpute","fcompress","regression")
print(tbl.true)

res[[exp.id]] = list()
res[[exp.id]]$tbl = tbl.true
res[[exp.id]]$fimpute = model.fimp
res[[exp.id]]$fpca = model.fpca
res[[exp.id]]$freg = model.fimp

# PLOT EXAMPLES
par(mfrow=c(2,2), cex=1.3)
ind = 1:2 + 10
idx = as.numeric(model.freg$id)

lims = c(min(ftrue[idx,][ind,],simulation$fobs[ind,],na.rm = TRUE) - 2,
         max(ftrue[idx,][ind,],simulation$fobs[ind,],na.rm = TRUE) + 2)
matplot(t(ftrue[idx,][ind,]),t='l',lty=1,lwd=4,ylim = lims,xlab = "time")
matplot(t(simulation$fobs[ind,]),t='p',lty=2,add=T,lwd=2,pch="x")
matplot(t(model.mean$fit[ind,]),t='l',lty=2,add=T,lwd=2)
title("Mean prediction")

matplot(t(ftrue[idx,][ind,]),t='l',lty=1,lwd=4,ylim = lims,xlab = "time")
matplot(t(simulation$fobs[ind,]),t='p',lty=2,add=T,lwd=2,pch="x")
matplot(t(model.fimp$fit[ind,]),t='l',lty=2,add=T,lwd=2)
title("Functional impute")

matplot(t(ftrue[idx,][ind,]),t='l',lty=1,lwd=4,ylim = lims,xlab = "time")
matplot(t(simulation$fobs[ind,]),t='p',lty=2,add=T,lwd=2,pch="x")
matplot(t(model.fpca$fit[ind,]),t='l',lty=2, add=T,lwd=2)
title("Functional PCA")

matplot(t(ftrue[idx,][ind,]),t='l',lty=1,lwd=4,ylim = lims,xlab = "time")
matplot(t(simulation$fobs[ind,]),t='p',lty=2,add=T,lwd=2,pch="x")
matplot(t(model.freg$fit[ind,]),t='l',lty=2, add=T,lwd=2)
title("Functional regression")

par(mfrow=c(1,2), cex=1.3)
matplot(t(-model.fimp$v)[,1:3],t='l',lwd = 4)
title("First 3 PCs from fPCA method")
matplot(t(model.fpca$v)[,1:3],t='l', lwd = 4)
title("First 3 singular vectors from SFI method")

par(mfrow=c(1,2), cex=1.3)
plot(cv.err ~ lambda, model.fimp$meta, type='o', ylab="CV error")
title("Sample error by cross-validation")
model.fimp$meta$cv.K = round(model.fimp$meta$cv.K)
plot(cv.K ~ lambda, model.fimp$meta, ylab="CV rank")
title("Rank by cross-validation")

n = nrow(model.fcmp$u)
#plot(model.fcmp$u[,1:2], col = 1 + (1:n > n/3), xlab = "Component 1", ylab="component 2")
#title("Decomposition of (Y(t),X_1(t),X_2(t))")
}

joint.tbl = res[[1]]$tbl
for (i in 2:length(res)){
  joint.tbl = joint.tbl + res[[i]]$tbl
}
joint.tbl = joint.tbl/length(res)
res[[1]]$fimpute$meta


# PLOT model.fimp
par(mfrow=c(1,1))
plot(model.fimp$meta[1:100,c(1,2)],ylim=c(0,10))
lines(model.fimp$meta[1:100,c(1,3)] )

joint.tbl
#cbPalette <- c("#8c1515", "#007c92")
cbPalette = list()
cbPalette[[1]] <- c("#8c1515", "#f4f4f4")
cbPalette[[2]] <- c("#f4f4f4", "#007c92")

stdPalette =  c("#8c1515", "#0098db","#eaab00", "#009b76", "#e98300", "#53284f", "#d2c295")

pid = 2

ids = c(1,3)
dd = data
dd$id = as.factor(dd$id)
dd$time = 5 + data$time*10
gg = as.numeric(colnames(model.fimp$Y))*10 + 5

pp = ggplot(aes(x = time, y = Y, color = id), data = dd[data$id %in% ids,]) +
  scale_fill_manual(values = cbPalette[[pid]]) + scale_colour_manual(values = cbPalette[[pid]]) +
  ylab("1st component") + xlab("age") +
  theme_set(theme_grey(base_size = 26)) + theme(legend.position="none", panel.background = element_rect(fill = "white",linetype = 1,colour = "grey50",size = 1,)) +
  scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0)) +
  stat_function(fun = approxfun(gg, simulation$basis[,1]) , color = stdPalette[1], size=1) +
  stat_function(fun = approxfun(gg, simulation$basis[,2]) , color = stdPalette[2], size=1) +
  stat_function(fun = approxfun(gg, simulation$basis[,3]) , color = stdPalette[3], size=1) +
  stat_function(fun = approxfun(gg, simulation$basis[,4]) , color = stdPalette[4], size=1) +
  stat_function(fun = approxfun(gg, simulation$basis[,5]) , color = stdPalette[5], size=1) +
  stat_function(fun = approxfun(gg, simulation$basis[,6]) , color = stdPalette[6], size=1) +
  stat_function(fun = approxfun(gg, simulation$basis[,7]) , color = stdPalette[7], size=1)
pp
ggsave(paste0("~/Dropbox/Presentations/Mobilize17/images/fcomplete/splines.pdf"))

pp = ggplot(aes(x = time, y = Y, color = id), data = dd[data$id %in% ids,]) +
  scale_fill_manual(values = cbPalette[[pid]]) + scale_colour_manual(values = cbPalette[[pid]]) +
  ylab("1st component") + xlab("age") +
  theme_set(theme_grey(base_size = 26)) + theme(legend.position="none", panel.background = element_rect(fill = "white",linetype = 1,colour = "grey50",size = 1,)) +
  scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0)) +
  stat_function(fun = approxfun(gg, model.fimp$v[1,]) , color = stdPalette[3], size=1) +
  stat_function(fun = approxfun(gg, model.fimp$v[2,]) , color = stdPalette[4], size=1) #+
  # stat_function(fun = approxfun(gg, simulation$basis[,3]) , color = stdPalette[3], size=1) +
  # stat_function(fun = approxfun(gg, simulation$basis[,4]) , color = stdPalette[4], size=1) +
  # stat_function(fun = approxfun(gg, simulation$basis[,5]) , color = stdPalette[5], size=1) +
  # stat_function(fun = approxfun(gg, simulation$basis[,6]) , color = stdPalette[6], size=1) +
  # stat_function(fun = approxfun(gg, simulation$basis[,7]) , color = stdPalette[7], size=1)
pp
ggsave(paste0("~/Dropbox/Presentations/Mobilize17/images/fcomplete/fpca-basis.pdf"))




pp = ggplot(aes(x = time, y = Y, color = id), data = dd[data$id %in% ids,]) +
  scale_fill_manual(values = cbPalette[[pid]]) + scale_colour_manual(values = cbPalette[[pid]]) +
  ylab("1st component") + xlab("age") +
  geom_point(size = 3) + theme_set(theme_grey(base_size = 26)) + theme(legend.position="none", panel.background = element_rect(fill = "white",linetype = 1,colour = "grey50",size = 1,)) +
  scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0)) +
  stat_function(fun = approxfun(lowess(dd$time,dd$Y)), size = 1, colour = "#555555") + xlim(5,15) + ylim(-4,4)

pp = pp +
  stat_function(fun = approxfun(gg, simulation$ftrue[ids[1],]) , color = cbPalette[[pid]][1], size=1, linetype='dashed') +
  stat_function(fun = approxfun(gg, simulation$ftrue[ids[2],]) , color = cbPalette[[pid]][2], size=1, linetype='dashed')
  # stat_function(fun = approxfun(gg, simulation$ftrue[ids[3],]) , color = cbPalette[3], size=1, linetype='dashed')
  # stat_function(fun = approxfun(gg, simulation$ftrue[4,]) , color = cbPalette[4], size=1, linetype='dashed')
pp
ggsave(paste0("~/Dropbox/Presentations/Mobilize17/images/fcomplete/observed-",pid,".pdf"))


pp + ggtitle("Individual mean") +
  stat_function(fun = approxfun(gg, model.mean$fit[ids[1],]) , color = cbPalette[[pid]][1], size=1.5, alpha = min(1,1.005 - (pid==2) )) +
  stat_function(fun = approxfun(gg, model.mean$fit[ids[2],]) , color = cbPalette[[pid]][2], size=1.5, alpha = min(1,1.005 - (pid==1) ))
  # stat_function(fun = approxfun(gg, model.mean$fit[ids[3],]) , color = cbPalette[3], size=1.5)
  # stat_function(fun = approxfun(gg, model.mean$fit[4,]) , color = cbPalette[4], size=1.5)
ggsave(paste0("~/Dropbox/Presentations/Mobilize17/images/fcomplete/2-curves-mean-",pid,".pdf"))

pp + ggtitle("Sparse PCA") +
  stat_function(fun = approxfun(gg, model.fpca$fit[ids[1],]) , color = cbPalette[[pid]][1], size=1.5, alpha = min(1,1.005 - (pid==2) )) +
  stat_function(fun = approxfun(gg, model.fpca$fit[ids[2],]) , color = cbPalette[[pid]][2], size=1.5, alpha = min(1,1.005 - (pid==1) ))
  # stat_function(fun = approxfun(gg, model.fimp$fit[ids[3],]) , color = cbPalette[3], size=1.5)
  # stat_function(fun = approxfun(gg, model.fimp$fit[4,]) , color = cbPalette[4], size=1.5)
ggsave(paste0("~/Dropbox/Presentations/Mobilize17/images/fcomplete/2-curves-fpca-",pid,".pdf"))

pp + ggtitle("Sparse Impute") +
  stat_function(fun = approxfun(gg, model.fimp$fit[ids[1],]) , color = cbPalette[[pid]][1], size=1.5, alpha = min(1,1.005 - (pid==2) )) +
  stat_function(fun = approxfun(gg, model.fimp$fit[ids[2],]) , color = cbPalette[[pid]][2], size=1.5, alpha = min(1,1.005 - (pid==1) ) )
  # stat_function(fun = approxfun(gg, model.fpca$fit[ids[3],]) , color = cbPalette[3], size=1.5)
  # stat_function(fun = approxfun(gg, model.fpca$fit[4,]) , color = cbPalette[4], size=1.5)
ggsave(paste0("~/Dropbox/Presentations/Mobilize17/images/fcomplete/2-curves-fimp-",pid,".pdf"))


