# Regenerate and reinstall fimpute
library("roxygen2") ; roxygenize()
library("devtools")
#devtools::install(".")
library("fcomplete")
library("ggplot2")
library("latex2exp")
library("parallel")

res = list()
dgrid = 51
d = 7

ns = c(100,500,1000)

# SIMULATE DATA
test.experiment = function(exp.id){
  res = list()
  for (i in ns){
  #set.seed(35 + exp.id)
  #simulation = fsimulate(dgrid = dgrid,clear = 0.9, n = 50, noise.mag = 0.3, d = d, K = 2)
  simulation = fsimulate(dgrid = dgrid,num_points = 3, n = i, noise.mag = 0.1, d = d, K = 4)
  data = simulation$data
  ftrue = simulation$ftrue
  K = simulation$params$K

  # TUNING PARAMS
  lambdas.pca = seq(0,0.1,length.out = 5)
  #lambdas.pca = 0
  lr = 0.5
  #lambdas.reg = seq(0,0.2,length.out = 20)
  #lambdas.reg = 0
  #lambdas.pca = 1

  model.mean = fregression(Y ~ time | id, data, method = "mean", bins = dgrid)

  ptm <- proc.time()
  model.fpca = fregression(Y ~ time | id, data, lambda = 0, K = 2:d, thresh = 1e-4, method = "fpcs", bins = dgrid, maxIter = 1000)
  time.fpca = (proc.time() - ptm)[3]

  ptm <- proc.time()
  model.fimp = fregression(Y ~ time | id, data, method="proximal_grad_grid", lambda = lambdas.pca, lr = lr, thresh = 1e-4, final = "soft", maxIter = 10000, fold = 5, cv.ratio = 0.05, bins = dgrid, verbose = 0)
  #model.fimp = fregression(Y ~ time | id, data, lambda = lambdas.pca, thresh = 1e-4, final = "soft", maxIter = 10000, fold = 5, cv.ratio = 0.05, bins = dgrid, verbose = 0)
  time.fimp = (proc.time() - ptm)[3]

  ptm <- proc.time()
  model.fimp.pg = fregression(Y ~ time | id, data, method="proximal_grad", lambda = lambdas.pca, lr = lr, thresh = 1e-4, final = "soft", maxIter = 10000, fold = 5, cv.ratio = 0.05, bins = dgrid, verbose = 0)
  time.fimp.pg = (proc.time() - ptm)[3]

  plot(lambdas.pca, model.fimp$loss)
  mean(((ftrue - model.fimp$fit)[!is.na(simulation$fobs)])**2)
  mean(((ftrue - model.fimp.pg$fit)[!is.na(simulation$fobs)])**2)

  #model.X1 = fregression(X1 + X2 + Y ~ time | id, data, K=4, lambda = lambdas.pca, thresh = 1e-10, method = "fimpute", bins = dgrid, maxIter = 5000)
  #model.fslr = fregression(Y ~ time + U1 + U2 + U3 + U4 | id, data, model.X1$u, K=3, thresh = 1e-10, lambda = lambdas.reg, method = "fimpute", bins = dgrid, maxIter = 5000)

  # REPORT RESULTS
  vv = mean((ftrue - mean(data$Y))**2)
  errors = c(
    mean((ftrue - model.fpca$fit)**2) / vv,
    mean((ftrue - model.fimp$fit)**2) / vv,
    mean((ftrue - model.fimp.pg$fit)**2) / vv
    #  mean((ftrue - model.fslr$fit)**2)
  )
  res[[i]] = list(errors = errors, time= c(time.fpca, time.fimp, time.fimp.pg))
  }
  res
}
res = mclapply(1:32, test.experiment, mc.cores = 8)
#save(res, file="simulation.study.Rda")

errors = c()
for (i in 1:length(res)){
  for (n in ns){
    errors = rbind(errors, c(n, res[[i]][[n]]$errors))
  }
}
times = c()
for (i in 1:length(res)){
  for (n in ns){
    times = rbind(times, c(n, res[[i]][[n]]$time))
  }
}
print(errors)
colMeans(errors)

colnames(errors) = c("n",
                     "fPCA",
                     "SLI",
                     "PG")
colnames(times) = colnames(errors)

library(tidyr)
errors.long = data.frame(errors) %>% gather(method, value, fPCA:PG)
errors.long = errors.long[errors.long$n==100,]

theme_set(theme_grey(base_size = 26))
dev.new()
p <- ggplot(errors.long, aes(x=method, y=value, color=method)) +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none", panel.background = element_rect(fill = "white",linetype = 1,colour = "grey50",size = 1,)) +
  ylab("NMSE") + xlab("") + ggtitle("Model performance") + ylim(c(0.3,0.7)) +
  geom_boxplot(size=1)
p
ggsave("figures/simulation-errors.pdf",width=7,height=5)

times.long = data.frame(times) %>% gather(method, value, fPCA:PG)
times.long = times.long[times.long$n==100,]
times.long$value = log10(times.long$value)
p <- ggplot(times.long, aes(x=method, y=value, color=method)) +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none", panel.background = element_rect(fill = "white",linetype = 1,colour = "grey50",size = 1,)) +
  ylab("time [log10(s)]") + xlab("") + ggtitle("Computation time") +
  geom_boxplot(size=1)
p
ggsave("figures/simulation-times.pdf",width=7,height=5)

library(dplyr)

times.long = data.frame(times) %>% gather(method, value, fPCA:PG)
times.long = times.long %>%
  group_by(n,method) %>%
  summarise(value = mean(value))
times.long$value = log10(times.long$value)

p <- ggplot(times.long, aes(x=n, y=value, color=method, background="white")) +
  theme(plot.title = element_text(hjust = 0.5), legend.key = element_rect(fill = "white"), panel.background = element_rect(fill = "white",linetype = 1,colour = "grey50",size = 1,)) +
  ylab("time [log10(s)]") + xlab("") + ggtitle("Computation time") +
  geom_line(size=1) + geom_point(size=1.5)
p
ggsave("figures/simulation-times-n.pdf",width=9,height=5)
