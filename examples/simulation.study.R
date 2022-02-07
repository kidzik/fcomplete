library("fcomplete")
library("ggplot2")
library("latex2exp")
library("parallel")
library("tidyr")
library("dplyr")

res = list()

dgrid = 51 # number of grid points
d = 7 # number of splines
ns = c(100) # number of subjects
nexperiments = 5 # set to 1000 for reproducing results from the paper
ncores = 5 # number of cores to use for the experiment

# Run one series of experiments with different number of subjects as in the ns variable
test.experiment = function(exp.id){
  set.seed(35 + exp.id)
  res = list()

  # iterate through different # of subject settings
  for (i in ns){
    # simulate data
    simulation = fsimulate(dgrid = dgrid,num_points = 3, n = i, noise.mag = 0.1, d = d, K = 4)
    data = simulation$data
    ftrue = simulation$ftrue
    K = simulation$params$K

    # Set tuning parameters, constant for all simulations
    lambdas.pca = seq(0,0.1,length.out = 5)
    lr = 0.5

    model.mean = fregression(Y ~ time | id, data, method = "mean", bins = dgrid)

    ptm <- proc.time()
    model.fpca = fregression(Y ~ time | id, data, lambda = 0, K = 2:d, thresh = 1e-4, method = "fpcs", bins = dgrid, maxIter = 1000)
    time.fpca = (proc.time() - ptm)[3]

    ptm <- proc.time()
    model.fimp = fregression(Y ~ time | id, data, method="pg_grid", lambda = lambdas.pca, lr = lr, thresh = 1e-4, final = "soft", maxIter = 10000, fold = 5, cv.ratio = 0.05, bins = dgrid, verbose = 0)
    time.fimp = (proc.time() - ptm)[3]

    ptm <- proc.time()
    model.fimp.pg = fregression(Y ~ time | id, data, method="pg", lambda = lambdas.pca, lr = lr, thresh = 1e-4, final = "soft", maxIter = 10000, fold = 1, cv.ratio = 0.05, bins = dgrid, verbose = 0)
    time.fimp.pg = (proc.time() - ptm)[3]

    # Collect results
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

# Run experiments
res = mclapply(1:nexperiments, test.experiment, mc.cores = ncores)

# Combine errors into one matrix
errors = c()
for (i in 1:length(res)){
  for (n in ns){
    errors = rbind(errors, c(n, res[[i]][[n]]$errors))
  }
}

# Combine computation time into one matrix
times = c()
for (i in 1:length(res)){
  for (n in ns){
    times = rbind(times, c(n, res[[i]][[n]]$time))
  }
}

# Display results
print(errors)

colnames(errors) = c("n", "fPCA", "SLI", "PG")
colnames(times) = colnames(errors)

errors.long = data.frame(errors) %>% gather(method, value, fPCA:PG)
errors.long = errors.long[errors.long$n==100,]

# Plot the performance figure
theme_set(theme_grey(base_size = 26))
dev.new()
p <- ggplot(errors.long, aes(x=method, y=value, color=method)) +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none", panel.background = element_rect(fill = "white",linetype = 1,colour = "grey50",size = 1,)) +
  ylab("NMSE") + xlab("") + ggtitle("Model performance") + ylim(c(0.3,0.7)) +
  geom_boxplot(size=1)
p
ggsave("figures/simulation-errors.pdf",width=7,height=5)

# Plot the computation time figure
times.long = data.frame(times) %>% gather(method, value, fPCA:PG)
times.long = times.long[times.long$n==100,]
times.long$value = log10(times.long$value)
p <- ggplot(times.long, aes(x=method, y=value, color=method)) +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none", panel.background = element_rect(fill = "white",linetype = 1,colour = "grey50",size = 1,)) +
  ylab("time [log10(s)]") + xlab("") + ggtitle("Computation time") +
  geom_boxplot(size=1)
p
ggsave("figures/simulation-times.pdf",width=7,height=5)

# Plot computation time as a function of n
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

