roxygen2::roxygenize(".")
#devtools::install(".")
library("fcomplete")
library("ggplot2")

spnbmd = ElemStatLearn::bone

model.proxgrad = fregression(as.formula(paste0("spnbmd ~ age | idnum")), spnbmd,
                             lambda= c(1,2), thresh = 1e-10, maxIter = 100,
                             method = "proximal_grad", final = "soft",
                             K=2, d=7, fold = 5, lr = 0.1)

model.impute = fregression(as.formula(paste0("+ spnbmd ~ age | idnum")), spnbmd,
                           lambda= 0.1, thresh = 1e-10, maxIter = 10000,
                           method = "fimpute", final = "soft",
                           K=2, d=7, fold = 5)

plot(model.impute, rows=3:5)
plot(model.proxgrad, rows=3:5)
fcomplete:::summary.fcomplete(model.impute)
#fcomplete:::predict.fcomplete(model.impute,ids=c(3,5),time=c(10,15,22))

## PCA
spnbmd$rnd = spnbmd$spnbmd + rnorm(length(spnbmd$spnbmd))*0.1
model.proxgrad = fregression(as.formula(paste0("spnbmd + rnd ~ age | idnum")), spnbmd,
                             lambda= c(1,2), thresh = 1e-10, maxIter = 100,
                             method = "proximal_grad", final = "soft",
                             K=2, d=7, fold = 5, lr = 0.1)


model.impute = fregression(as.formula(paste0("spnbmd + rnd ~ age | idnum")), spnbmd,
                             lambda= c(1,2), thresh = 1e-10, maxIter = 100,
                             final = "soft",
                             K=2, d=7, fold = 5, lr = 0.1)



predict(model.impute, newdata=spnbmd)
coef(model.impute)
residuals(model.impute)
fitted(model.impute)

predict(model.impute, ids=c(2,5))
predict(model.impute, ids=c(2,5), time= c(10,20))

ind = 50 + 1:3
#plt = plot_preds(model.impute$Y[ind,], NULL, model.impute$fit[ind,],
#                 filename="pred-freg-data.pdf", title = "Sparse Functional Impute")

model.fpca = fregression(as.formula(paste0("spnbmd ~ age | idnum")), spnbmd,
                           method = "fpcs", d=7, K=2)
predict(model.fpca, newdata=spnbmd)
coef(model.fpca)
residuals(model.fpca)
fitted.values(model.fpca)

ind = 50 + 1:3
#plt = plot_preds(model.fpca$Y[ind,], NULL, model.fpca$fit[ind,],
#                 filename="pred-freg-data.pdf", title = "Sparse Functional Impute")
