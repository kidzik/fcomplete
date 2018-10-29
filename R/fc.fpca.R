# TODO: This file is ugly for now
# The main purpose is to wrap fpcs functions to our framework
# @export
fc.fpca = function(X, d=5, K=2, grid.l=0:99/99){
  ids = X[,1]
  time = as.numeric(X[,2])
  values = as.numeric(X[,3])

  # X[is.na(X[,3]),3] = mean(X[,3],na.rm = TRUE)

  mint = min(time)
  maxt = max(time)

  mapto01 = function(time) { (time - mint) / (maxt - mint) }
  mapfrom01 = function(time) {  mint + time * (maxt - mint) }

  time = mapto01(time)

  ## candidate models for fitting
  #  M.set<-c(4,5,6)
  M.set<-d
  #  r.set<-c(2,3,4)
  r.set<-K
  ##parameters for fpca.mle
  ini.method="EM"
  basis.method="bs"
  sl.v=rep(0.5,10)
  max.step=50
  #  grid.l= min(X[,2]) + (max(X[,2]) - min(X[,2])) * 0:99/99
  grids = (grid.l - (min(grid.l))) / (max(grid.l) - min(grid.l)) #seq(0,1,0.002)
  ##fit candidate models by fpca.mle
  result=fpca.mle(cbind(ids,values,time), M.set,r.set,ini.method, basis.method,sl.v,max.step,grid.l,grids)
  muest<-result$fitted_mean

  #get predicted trajectories on a fine grid: the same grid for which mean and eigenfunctions are evaluated
  ##rescaled grid
  grids.new<-result$grid
  M<-result$selected_model[1]
  r<-result$selected_model[2]
  evalest<-result$eigenvalues     ## estimated
  sig2est<-result$error_var        ## estimated
  eigenfest<-result$eigenfunctions

  #X[,c(1,3,2)]
  fpcs<-fpca.score.fixed(cbind(ids,values,time),grids.new,muest,evalest,eigenfest,sig2est,r)
  list(fit = t(fpca.pred(fpcs, muest, eigenfest)), sigma.est = sqrt(sig2est), mu.est = muest, selected_model = result$selected_model, fpcs = fpcs, v = eigenfest)
}
