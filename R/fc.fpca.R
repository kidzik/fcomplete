#########
fpca.score.fixed<-function(data.m,grids.u,muhat,eigenvals,eigenfuncs,sig2hat,K){
  ##estimated conditional principal component scores (BLUPs): \hat{E(\xi_k|Y)}, k=1,...,K
  ##Name:FPcScore
  ##para:
  ##     data.m -- data matrix; same as input for fpca.mle
  ##     grids.u -- grid of time points used in evaluating the mean and eigenfunctions (on the original scale); (returned by fpca. mle)
  ##     muhat,eigenvals, eigenfuncs, sig2hat -- (estimated) mean, eigenvalues, eigenfunctions and noise variance; (returned by fpca.mle)
  ##     K -- number of eigenfunctions used in the model, i.e., (estimated) dimension of the process
  ##return: first K conditional PC scores (the BLUP estimates): n by K
  temp<-table(as.numeric(data.m[,1]))
  n<-length(temp)             ##     number of curves;
  m.l<-as.vector(temp)        ##     m.l -- number of time points per curve
  result<-matrix(0,n,K)       ##First K FPC scores for each subject

  N <- length(grids.u)        ## number of time points on the grid
  evalmat <- diag(eigenvals[1:K])  ## diagonal matrix of the first K (estimated) eigenvalues
  current <- 0  ## current index
  eigenfuncs.u<-t(eigenfuncs)   ## dimmension: grid_length by K

  data.u<-matrix(as.numeric(as.vector(as.matrix(data.m[,-1]))),nrow=nrow(data.m[,-1]),ncol=ncol(data.m[,-1]))     ##convert obs matrix to be numierc

  for (i in 1:n){
    Y <- as.vector(data.u[(current+1):(current+m.l[i]),1])  ## observed  measurements of ith curve
    meastime <- data.u[(current+1):(current+m.l[i]),2] ## measurement times of the ith curve
    gridtime <- ceiling(N*meastime)   ## project measurement time onto the grid
    gridtime[gridtime <= 0] = 1 #quickfix
    muy <- muhat[gridtime]
    Phiy  <- matrix(eigenfuncs.u[gridtime,1:K],ncol=K)
    Sigy <- Phiy %*% evalmat %*% t(Phiy) + sig2hat * diag(m.l[i])
    temp.y<-matrix(Y-muy)
    result[i,] <- evalmat %*% t(Phiy) %*% solve(Sigy,temp.y)
    current <- current + m.l[i]
  }
  row.names(result) = levels(data.m[,1])[as.numeric(names(temp))]
  result
#  list(scores = result, ids = ] )
}


#' @export
fc.fpca = function(X, d=5, K=2, grid.l=0:99/99){
  X[is.na(X[,3]),3] = mean(X[,3],na.rm = TRUE)

  mint = min(X[,2])
  maxt = max(X[,2])

  mapto01 = function(time) { (time - mint) / (maxt - mint) }
  mapfrom01 = function(time) {  mint + time * (maxt - mint) }

  X[,2] = mapto01(X[,2])

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
  result=fpca.mle(X[,c(1,3,2)], M.set,r.set,ini.method, basis.method,sl.v,max.step,grid.l,grids)
  muest<-result$fitted_mean

  #get predicted trajectories on a fine grid: the same grid for which mean and eigenfunctions are evaluated
  ##rescaled grid
  grids.new<-result$grid
  M<-result$selected_model[1]
  r<-result$selected_model[2]
  evalest<-result$eigenvalues     ## estimated
  sig2est<-result$error_var        ## estimated
  eigenfest<-result$eigenfunctions

  fpcs<-fpca.score.fixed(X[,c(1,3,2)],grids.new,muest,evalest,eigenfest,sig2est,r)
  list(fit = t(fpca.pred(fpcs, muest, eigenfest)), sigma.est = sqrt(sig2est), mu.est = muest, selected_model = result$selected_model, fpcs = fpcs)
}
