library("fcomplete")
library("ggplot2")

spnbmd = ElemStatLearn::bone

model.impute = fregression(as.formula(paste0("spnbmd:age ~ 1 | idnum")), spnbmd,
                           lambda= 0:10/10, thresh = 1e-10, maxIter = 10000,
                           method = "fimpute", final = "soft",
                           K=2, d=7, fold = 5)

ind = 50 + 1:3
plt = plot_preds(model.impute$Y[ind,], NULL, model.impute$fit[ind,],
                 filename="pred-freg-data", title = "Sparse Functional Impute")

model.fpca = fregression(as.formula(paste0("spnbmd:age ~ 1 | idnum")), spnbmd,
                           method = "fpcs", d=7, K=2)

ind = 50 + 1:3
plt = plot_preds(model.fpca$Y[ind,], NULL, model.fpca$fit[ind,],
                 filename="pred-freg-data", title = "Sparse Functional Impute")
