library("fcomplete")
library("ggplot2")
#library("latex2exp")
library("parallel")
library("tidyr")
library("dplyr")
library("fdapace")
#library("refund")
library("face")
library("mfaces")

res = list()

dgrid = 51 # number of grid points
d = 7 # number of splines
ns = c(100, 500, 1000, 3000) # number of subjects
nexperiments = 10 # set to 1000 for reproducing results from the paper
ncores = 10 # number of cores to use for the experiment

# Load data
load("data/gdi.Rda")

sampling.data = all.data.filtered[,c("Patient_ID","age")]
sampling.data$age = (sampling.data$age - min(sampling.data$age)) / (max(sampling.data$age) - min(sampling.data$age))

data.sampling = function(){
  sids = table(sampling.data$Patient_ID)
  sid = sample(as.numeric(names(sids[sids>=4])))[1]
  sort(sampling.data$age[sampling.data$Patient_ID == sid])
}

resample.data = function(simulation){
  simulation$fobs[] = NA
  n = nrow(simulation$fobs)
  T = ncol(simulation$fobs)

  ids = c()
  vals = c()
  time = c()

  for (i in 1:n){
    timepoints = data.sampling()

    for (tp in timepoints){
      ids = c(ids, i)
      k = which(1:T/T >= tp)[1]
      tt = k/T
      time = c(time, tt)
      vals = c(vals, simulation$fnoisy[i,k])
      simulation$fobs[i,k] = simulation$fnoisy[i,k]
    }
  }
  simulation$data = data.frame(id = ids, time = time, Y = vals)
  simulation
}


# Run one series of experiments with different number of subjects as in the ns variable
test.experiment = function(exp.id){
  set.seed(30 + exp.id)
  res = list()

  # iterate through different # of subject settings
  for (i in ns){
    # simulate data
    simulation = fsimulate(dgrid = dgrid,num_points = 3, n = i, noise.mag = 0.1, d = d, K = 4)
#    simulation = resample.data(simulation)

    data = simulation$data
    ftrue = simulation$ftrue
    K = simulation$params$K

    # Set tuning parameters, constant for all simulations
    lambdas.pca = seq(0,0.1,length.out = 5)
    lr = 0.5

    model.mean = fregression(Y ~ time | id, data, method = "mean", bins = dgrid)

    ptm <- proc.time()
    model.fimp.pg.multi = fregression(Y + X1 + X2 ~ time | id, data, method="fimpute", lambda = lambdas.pca, lr = lr, thresh = 1e-4, final = "soft", maxIter = 10000, fold = 5, cv.ratio = 0.05, bins = dgrid, verbose = 0)
    time.fimp.mv = (proc.time() - ptm)[3]

    ptm <- proc.time()
    model.fpca = fregression(Y ~ time | id, data, lambda = 0, K = 2:d, thresh = 1e-4, method = "fpcs", bins = dgrid, maxIter = 1000)
    time.fpca = (proc.time() - ptm)[3]

    ptm <- proc.time()
    model.fimp = fregression(Y ~ time | id, data, method="fimpute", lambda = lambdas.pca, lr = lr, thresh = 1e-4, final = "soft", maxIter = 10000, fold = 5, cv.ratio = 0.05, bins = dgrid, verbose = 0)
    time.fimp = (proc.time() - ptm)[3]

    ptm <- proc.time()
    model.fimp.pg = fregression(Y ~ time | id, data, method="pg_grid", lambda = lambdas.pca, lr = lr, thresh = 1e-4, final = "soft", maxIter = 10000, fold = 1, cv.ratio = 0.05, bins = dgrid, verbose = 0)
    time.fimp.pg = (proc.time() - ptm)[3]

    ptm <- proc.time()
    n = max(simulation$data$id)
    Ly = list()
    t = list()
    for (j in 1:n){
      p = simulation$data[simulation$data$id==j,]
      Ly[[j]]=p$Y
      t[[j]]=p$time
    }
    fit.fdapace = FPCA(Ly, t)
    time.fimp.fdapace = (proc.time() - ptm)[3]

    # face prepare data
    face.data = list()
    for (j in 1:3){
      face.data[[j]] = data[,c(1,2,2+j)]
      colnames(face.data[[j]]) = c("subj","argvals","y")
    }

    n = max(face.data[[1]]$subj)
    grid = rep(1,n) %*% t(c(simulation$params$grid))
    newdata = cbind(rep(1:n, each = 51), c(t(grid)))

    newdata.df = list()

    for (j in 1:3){
      newdata.df[[j]] = data.frame(subj = newdata[,1], argvals = newdata[,2])
      newdata.df[[j]]$y = NA
      newdata.df[[j]] = rbind(newdata.df[[j]], face.data[[j]])
    }

    # face single
    ptm <- proc.time()
    face.fit = face.sparse(face.data[[1]])
    m = dim(grid)[2]
    preds = predict(face.fit, newdata.df[[1]])
    fit.face = t(matrix(preds$y.pred[1:(n*m)], m))

    time.fimp.face = (proc.time() - ptm)[3]

    # face multi
    ptm <- proc.time()
    face.fit.multi = mfaces::mface.sparse(face.data)
    preds = predict(face.fit.multi, newdata.df)
    time.fimp.face.mv = (proc.time() - ptm)[3]

    fit.face.mv = t(matrix(preds$y.pred[[1]][1:(n*m)], m))

    # Collect results
    vv = mean((ftrue - mean(data$Y))**2)

    errors = c(
      mean((ftrue - model.fpca$fit)**2) / vv,
      mean((ftrue - model.fimp$fit)**2) / vv,
      mean((ftrue - model.fimp.pg$fit)**2) / vv,
      mean((ftrue - fitted(fit.fdapace))**2) / vv,
      mean((ftrue - fit.face)**2) / vv,
      mean((ftrue - model.fimp.pg.multi$multiFit[[1]])**2)/vv,
      mean((ftrue - fit.face.mv)**2) / vv
    )
    time = c(time.fpca,
             time.fimp,
             time.fimp.pg,
             time.fimp.fdapace,
             time.fimp.face,
             time.fimp.mv,
             time.fimp.face.mv
    )
    res[[i]] = list(errors = errors, time= time)
    res[[i]]
  }
  res
}

# Run experiments
res = mclapply(1:nexperiments, test.experiment, mc.cores = ncores)
save(res, file="experiment-3.Rda")

# Combine errors into one matrix
errors = c()
for (i in 1:length(res)){
  for (n in ns){
    if (typeof(res[[i]]) != "character")
      errors = rbind(errors, c(n, res[[i]][[n]]$errors))
  }
}

# Combine computation time into one matrix
times = c()
for (i in 1:length(res)){
  for (n in ns){
    if (typeof(res[[i]]) != "character")
      times = rbind(times, c(n, res[[i]][[n]]$time))
  }
}

# Display results
colMeans(errors)

colnames(errors) = c("n", "fPCA", "SLI", "PG", "fda.pace", "face", "SLR", "face.mv")
colnames(times) = colnames(errors)

errors.long = data.frame(errors) %>% gather(method, value, fPCA:face)
#errors.long = data.frame(errors) %>% gather(method, value, SLR:face.mv)
errors.long = errors.long[errors.long$n==100,]

# Plot the performance figure
theme_set(theme_classic(base_size = 26))
p <- ggplot(errors.long, aes(x=method, y=value, color=method)) +
  theme(axis.line.x.bottom=element_line(size=0.5), axis.line.y.left=element_line(size=0.5),
        plot.title = element_text(hjust = 0.5), legend.position="none",
        panel.background = element_rect(fill = "white",linetype = 0,colour = "grey50",size = 1,)) +
  ylab("NMSE") + xlab("") + ggtitle("Model performance") + ylim(c(min(errors.long$value),max(errors.long$value))) +
  geom_boxplot(size=1)
p
ggsave("figures/simulation-errors.pdf",width=7,height=5)

# Plot the computation time figure
times.long = data.frame(times) %>% gather(method, value, fPCA:face.mv)
# times.long = data.frame(times) %>% gather(method, value, SLR:face.mv)
times.long = times.long[times.long$n==100,]
times.long$value = log10(times.long$value)
p <- ggplot(times.long, aes(x=method, y=value, color=method)) +
  theme(axis.line.x.bottom=element_line(size=0.5), axis.line.y.left=element_line(size=0.5),
        plot.title = element_text(hjust = 0.5), legend.position="none",
        panel.background = element_rect(fill = "white",linetype = 0,colour = "grey50",size = 1,)) +
  ylab("time [log10(s)]") + xlab("") + ggtitle("Computation time") +
  geom_boxplot(size=1)
p
ggsave("figures/simulation-times.pdf",width=7,height=5)

# Plot computation time as a function of n
times.long = data.frame(times) %>% gather(method, value, fPCA:face)
times.long = times.long %>%
  group_by(n,method) %>%
  summarise(value = mean(value))
times.long$value = log10(times.long$value)

p <- ggplot(times.long, aes(x=n, y=value, color=method, background="white")) +
  theme(axis.line.x.bottom=element_line(size=0.5), axis.line.y.left=element_line(size=0.5),
        plot.title = element_text(hjust = 0.5),
        legend.key = element_rect(fill = "white", linewidth = 0, linetype = 0),
        panel.background = element_rect(fill = "white",linetype = 0,colour = "grey50",size = 1,)) +
  ylab("time [log10(s)]") + xlab("") + ggtitle("Computation time") +
  geom_line(size=1) + geom_point(size=1.5)
p
ggsave("figures/simulation-times-n.pdf",width=9,height=5)



