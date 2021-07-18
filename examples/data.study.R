library("fcomplete")
library("ggplot2")
library("parallel")
library("metafolio")
library("ggplot2")
library("ggthemes")
library("RColorBrewer")

sample_ratio = 1 # 1 for results from the study

clean_theme = theme_classic() + theme(text = element_text(size=18), panel.border = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))

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
  trialInfo = read.csv("/home/kidzik/Dropbox/DATA/CP/trialInfo_CP.csv")
}

all.data.subset = all.data[all.data$side == "L",]
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

pats = sample(pats)[1:floor(length(pats)*sample_ratio)]
all.data.filtered.sample = all.data.filtered[all.data.filtered$Patient_ID %in% pats,]

#################
# RUN CROSS-VAL #
#################
experiment.data = function(i)
{
  set.seed(i+100)
  # Sample data for testing

  var = "GDI"
  data = sample.long(all.data.filtered.sample, "Patient_ID", "age", var, ratio = 0.05, min.per.sbj = 3)

  # Set up parameters
  lambdas = list()
  lambdas[["GDI"]] = seq(4,9,length.out = 6)
  lambdas[["bmi"]] = seq(1,2,length.out = 5)
  d = 6
  K = 2

  # IMPUTE
  model.impute = fregression(as.formula(paste0(var," ~ age | Patient_ID")), data$train,
                             lambda= lambdas[["GDI"]], thresh = 1e-4, maxIter = 10000,
                             method = "fimpute", final = "soft",
                             K=K, d=d, fold = 5)

  model.impute.fpcs = fregression(as.formula(paste0(var," ~ age | Patient_ID")), data$train,
                                  thresh = 1e-4, method = "fpcs", K=K, d=d)

  # model.covariates = fregression(GDI + O2cost + speed ~ age  | Patient_ID, data$train,
  #                                method = "fimpute", thresh=1e-10, maxIter = 5000,
  #                                lambda = lambdas[[var]] / sqrt(120),
  #                                d=d,
  #                                K = 2)

  # model.regression = fregression(GDI ~ age + U1 | Patient_ID, data$train, model.covariates$u,
  #                                method = "fimpute", thresh=1e-10, maxIter = 5000,
  #                                K = 2,
  #                                lambda = 1,
  #                                d=d)

  errors = c(#mean((model.regression$fit - data$test.matrix)**2, na.rm = TRUE),
             mean((model.impute$fit - data$test.matrix)**2, na.rm = TRUE),
             mean((model.impute.fpcs$fit - data$test.matrix)**2, na.rm = TRUE),
             mean((mean(data$train.matrix,na.rm=TRUE) - data$test.matrix)**2, na.rm = TRUE)
  )

  names(errors) = c(#"reg",
                    "SLI",
                    "fPCA",
                    "mean")
  print(errors)
  list(errors = errors,
       model.fimp = model.impute,
       model.fpca = model.impute.fpcs,
       #model.fslr = model.regression,
       #       model.mean = model.mean,
       data = data)
}
models = mclapply(1:20, experiment.data, mc.cores = 20)



plot(models[[1]]$model.fimp$fit[!is.na(models[[1]]$data$test.matrix)], models[[1]]$model.fimp$Y[!is.na(models[[1]]$data$test.matrix)])

pred = models[[1]]$model.fimp$u[] %*% models[[1]]$model.fimp$d[1] %*% models[[1]]$model.fimp$v[1,]
mean((pred - models[[1]]$data$test.matrix)**2, na.rm = TRUE)
models[[1]]$errors

#df.quadtmp = data.frame(Patient_id = all.data$Patient_ID, quadriplegic = all.data$isQuadriplegic, triplegic = all.data$isTriplegic)
#df.quad = aggregate(list(quadriplegic = df.quadtmp$quadriplegic, triplegic = df.quadtmp$triplegic), list(df.quadtmp$Patient_id), max)
df.quad = aggregate(list(dxmod = as.character(trialInfo$dxmod)), list(Patient_id = trialInfo$Patient_ID), function(x){unique(x)[1]})
#names(df.quad) = c("Patient_id", "quadriplegic_bin","triplegic_bin")
df = data.frame(Patient_id = row.names(models[[1]]$data$train.matrix), trend.pc1 = models[[1]]$model.fimp$u[,1])
df.quad$dxmod = as.character(df.quad$dxmod)
df.quad = df.quad[df.quad$dxmod != "Femoral anteversion",]
df.quad = df.quad[df.quad$dxmod != "Hemiplegia type I",]

df.merged = merge(df,df.quad)
summary(lm(df.merged$trend.pc1 ~ df.merged$dxmod))

plt <- ggplot(df.merged, aes(x=df.merged$dxmod, y=df.merged$trend.pc1)) +
  geom_boxplot() +
  clean_theme +
  xlab("Subtype of the paralysis") +
  ylab("Score of the first PC")
print(plt)
ggsave(plt, filename = "figures/paralysis-subtypes.pdf",width = 10, height = 6)

age.min = min(all.data.filtered.sample$age)
age.max = max(all.data.filtered.sample$age)
age = age.min + (age.max - age.min) * ((0:50)/50)

## Figure components in data study
df1 = data.frame(age = age, GDI = models[[1]]$model.fimp$v[1,], component = "1")
df2 = data.frame(age = age, GDI = models[[1]]$model.fimp$v[2,], component = "2")
df = rbind(df1,df2)
plt = ggplot(df, aes(x = age, y = GDI, group = component, color = component)) + clean_theme + geom_line(size=1.5)
print(plt)
ggsave(plt, filename = "figures/data-components.pdf",width = 6, height = 4)

## FPCA
df1 = data.frame(age = age, GDI = models[[1]]$model.fpca$v[1,], component = "1")
df2 = data.frame(age = age, GDI = models[[1]]$model.fpca$v[2,], component = "2")
df = rbind(df1,df2)
plt = ggplot(df, aes(x = age, y = GDI, group = component, color = component)) + clean_theme + geom_line(size=1.5)
print(plt)
ggsave(plt, filename = "figures/data-components-FPCA.pdf",width = 6, height = 4)

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

# Lambdas
models[[exp.id]]$model.fimp$meta
models[[exp.id]]$model.fpca$selected_model

#############################
# OTHER PLOTS FOR THE PAPER #
#############################

dd = all.data.filtered[,c("Patient_ID","age","bmi","GDI")]
dd$Patient_ID = as.factor(dd$Patient_ID)
dd = dd[dd$bmi > 10,]

data_highlight=dd[dd$Patient_ID,]

# Figure 1: BMI over time
theme_set(theme_grey(base_size = 26))
pp = ggplot(aes(x = age, y = GDI, color = Patient_ID), data = dd[1:200,]) + ylab("GDI") + theme(legend.position="none", panel.background = element_rect(fill = "white",linetype = 1,colour = "grey50",size = 1,)) +
  stat_function(fun = approxfun(lowess(dd$age,dd$GDI)), size = 1.5, colour = "#000000")+ scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0))
pp + geom_point(size = 1.5,alpha=0.5)
ggsave("figures/data-points.pdf",width=7,height=5)

pp1 = pp + geom_line(size=0.7,alpha=0.05) + geom_point(size = 1.5, alpha=0.05)
pats = c(3995, 4241, 4215)
colors = c("#ff0000","#00ff00","#0000ff")
for (i in 1:3){
  pp1 = pp1 + geom_line(size=1,alpha=0.9, data = dd[dd$Patient_ID == pats[i], ], colour = colors[i]) + geom_point(size = 1.5, alpha=0.9, data = dd[dd$Patient_ID == pats[i], ], colour = colors[i])
}
pp1
ggsave("figures/data-grouped.pdf",width=7,height=5)

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
