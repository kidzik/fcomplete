library("fcomplete")
library("ggplot2")
library("parallel")
#library("metafolio")
library("ggplot2")
library("ggthemes")
library("RColorBrewer")

sample_ratio = 1 # sample patients from the  full dataset: 1 for results from the study, smaller numbers for quick tests
nexperiments = 2 # number of resampling experiments (100 for results in the paper)
ncores = 2 # number of cores for parallel runs

# theme for ggplot
clean_theme = theme_classic() + theme(text = element_text(size=18), panel.border = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))

# Load data
load("data/gdi.Rda")

# Resample data if needed
pats = as.numeric(names(table(all.data.filtered$Patient_ID)))
pats = sample(pats)[1:floor(length(pats)*sample_ratio)]
all.data.filtered.sample = all.data.filtered[all.data.filtered$Patient_ID %in% pats,]

# Data experiment with resampling
experiment.data = function(i)
{
  set.seed(i+1000)

  # Sample data for testing
  var = "GDI"
  data = sample.long(all.data.filtered.sample, "Patient_ID", "age", var, ratio = 0.05, min.per.sbj = 3)

  # Set up parameters
  lambdas = seq(1,20,length.out = 10)
  d = 6
  K = 5

  # IMPUTE
  model.impute = fregression(as.formula(paste0(var," ~ age | Patient_ID")), data$train,
                             lambda = lambdas, thresh = 1e-4, maxIter = 10000,
                             method = "fimpute", final = "soft",
                             K=K, d=d, fold = 5)

  model.impute.fpcs = fregression(as.formula(paste0(var," ~ age | Patient_ID")), data$train,
                                  thresh = 1e-4, method = "fpcs", K=2, d=d)

  errors = c(mean((model.impute$fit - data$test.matrix)**2, na.rm = TRUE),
             mean((model.impute.fpcs$fit - data$test.matrix)**2, na.rm = TRUE),
             mean((mean(data$train.matrix,na.rm=TRUE) - data$test.matrix)**2, na.rm = TRUE)
  )

  names(errors) = c("SLI",
                    "fPCA",
                    "mean")
  print(errors)
  list(errors = errors,
       model.fimp = model.impute,
       model.fpca = model.impute.fpcs,
       data = data)
}
models = mclapply(1:nexperiments, experiment.data, mc.cores = ncores)

# Analyze relation between the first PC and the subtype of the disease
df.quad = aggregate(list(dxmod = as.character(impairement.info$dxmod)), list(Patient_id = impairement.info$Patient_ID), function(x){unique(x)[1]})
df = data.frame(Patient_id = row.names(models[[1]]$data$train.matrix), trend.pc1 = models[[1]]$model.fimp$u[,1])
df.quad$dxmod = as.character(df.quad$dxmod)
df.quad = df.quad[df.quad$dxmod != "Femoral anteversion",]
df.quad = df.quad[df.quad$dxmod != "Hemiplegia type I",]
df.merged = merge(df,df.quad)
summary(lm(df.merged$`trend.pc1` ~ df.merged$`dxmod`))

plt <- ggplot(df.merged, aes(x=dxmod, y=trend.pc1)) +
  geom_boxplot() +
  clean_theme +
  xlab("Subtype of the paralysis") +
  ylab("Score of the first PC")
print(plt)
ggsave(plt, filename = "figures/paralysis-subtypes.pdf",width = 10, height = 6)

age.min = min(all.data.filtered.sample$age)
age.max = max(all.data.filtered.sample$age)
age = age.min + (age.max - age.min) * ((0:50)/50)

# Plot the first two components
df1 = data.frame(age = age, GDI = models[[1]]$model.fimp$v[1,], component = "PC1")
df2 = data.frame(age = age, GDI = models[[1]]$model.fimp$v[2,], component = "PC2")
df = rbind(df1,df2)
plt = ggplot(df, aes(x = age, y = GDI, group = component, color = component)) + clean_theme + geom_line(size=1.5)
print(plt)
ggsave(plt, filename = "figures/data-components.pdf",width = 6, height = 4)

# Plot how components affect the mean
df1 = data.frame(age = age, GDI = models[[1]]$model.fimp$cmeans + 10*models[[1]]$model.fimp$v[1,], curve = "mean + PC1")
df2 = data.frame(age = age, GDI = models[[1]]$model.fimp$cmeans, curve = "mean")
df3 = data.frame(age = age, GDI = models[[1]]$model.fimp$cmeans - 10*models[[1]]$model.fimp$v[1,], curve = "mean - PC1")
df = rbind(df1,df2,df3)
colors = gg_color_hue(3)
plt = ggplot(df, aes(x = age, y = GDI, group = curve, color = curve)) + clean_theme +
  geom_line(aes(linetype=curve, color=curve), size=1.5) +
  scale_linetype_manual(values=c("solid", "dashed", "dashed")) +
  scale_color_manual(values=c("black", colors[1], colors[2]))
print(plt)
ggsave(plt, filename = "figures/data-components-added.pdf",width = 6, height = 4)

# Summarize results
res = c()
for (i in 1:length(models)){
  if (is.null(attr(models[[i]],"class")))
    res = cbind(res, models[[i]]$errors)
}
rownames(res) = c("SLI","fPCA","mean")
cbind(rowMeans(res),
apply(res,FUN=sd,1))
colnames(res) = paste("run",1:length(models))
ind = row.names(models[[i]]$data$test.matrix) %in% models[[i]]$data$X$Patient_ID[models[[i]]$data$test.ob[1:3]]

## Additional plots for the paper
dd = all.data.filtered[,c("Patient_ID","age","bmi","GDI")]
dd$Patient_ID = as.factor(dd$Patient_ID)
dd = dd[dd$bmi > 10,]

# GDI over time
theme_set(theme_classic(base_size = 26))
pp = ggplot(aes(x = age, y = GDI, color = Patient_ID), data = dd[1:200,]) + ylab("GDI") +
  theme(axis.line.x.bottom=element_line(size=0.5), axis.line.y.left=element_line(size=0.5),
        legend.position="none", panel.background = element_rect(fill = "white",linetype = 0,colour = "grey50",size = 1,)) +
  stat_function(fun = approxfun(lowess(dd$age,dd$GDI)), size = 1.5, colour = "#000000")+ scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0))
pp + geom_point(size = 1.5,alpha=0.5)
ggsave("figures/data-points.pdf",width=7,height=5)

# GDI over time groupped by patient id
pp1 = pp + geom_line(size=0.7,alpha=0.05) + geom_point(size = 1.5, alpha=0.05)
pats = c(3995, 4241, 4215)
colors = c("#ff0000","#00ff00","#0000ff")
for (i in 1:3){
  pp1 = pp1 + geom_line(size=1,alpha=0.9, data = dd[dd$Patient_ID == pats[i], ], colour = colors[i]) + geom_point(size = 1.5, alpha=0.9, data = dd[dd$Patient_ID == pats[i], ], colour = colors[i])
}
pp1
ggsave("figures/data-grouped.pdf",width=7,height=5)


# Review: Data generation process

Ds = c()
for (i in 1:nexperiments){
  Ds = rbind(Ds, models[[i]]$model.fimp$d)
}
cm = colMeans(Ds)

pdf("figures/revision-decomposition.pdf",height = 5, width = 8)
boxplot(Ds/cm[1])
title("Diagonal of the spectral decomposition",xlab = "index",ylab="normalized value")
svals = ((1/3) * c(1, 0.6, 0.3, 0.2 *exp(-3), 0.1 *exp(-4)) + (2/3) * c(1.3,0.4, 0.2, 0.2 *exp(-3), 0.1 *exp(-4)))
points(  svals/svals[1], col="blue", lwd=7)
lines(  svals/svals[1], col="blue")
dev.off()
