library("devtools")
install(".")
source("tests/plot.helpers.R")

# Load simulated 'toy' data
data = read.table("examples/simulated.Rda")

# Print the data
head(data)

# Run our model
model = fregression(Y:time ~ 1|id, data, bins = 41, K = 3)

# Make some diagnostics plots
idx = 1:5
plot_preds(model$Y[idx,], model$fit[idx,], NULL,
           filename="pred-mean", title = "True curves & observations")

# Some output elements
model$v # principal components (elements from K + 1 are meaningless, I'll remove them)
model$u # loadings on the first K principal components

# 'patients' plotted on their 2D progression plane (third pattern skipped)
plot(model$u[,1:2])

# 'Progression' trends
matplot(t(model$v[1:3,]),t='l')
