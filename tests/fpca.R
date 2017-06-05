############################## example I: "easy" case
##load data
data(easy)
##sample trajectory
plot(easy$data[,3],easy$data[,2],xlab="grid",ylab="",type='p')
for(i in 1:500){
  cur<-easy$data[easy$data[,1]==i,]
  points(cur[,3],cur[,2],type="l")
}
## candidate models for fitting
M.set<-c(4,5,6)
r.set<-c(2,3,4)
##parameters for fpca.mle
ini.method="EM"
basis.method="bs"
sl.v=rep(0.5,10)
max.step=50
grid.l=seq(0,1,0.01)
grids=seq(0,1,0.002)
##fit candidate models by fpca.mle
result<-fpca.mle(easy$data, M.set,r.set,ini.method, basis.method,sl.v,max.step,grid.l,grids)
summary(result)
##rescaled grid

grids.new<-result$grid
##model selection result: the true model M=5, r=3 is selected with the smallest CV score among all converged models
M<-result$selected_model[1]
r<-result$selected_model[2]
##compare estimated eigenvalues with the truth
evalest<-result$eigenvalues     ## estimated
easy$eigenvalues    ## true
##compare estimated error variance with the truth
sig2est<-result$error_var        ## estimated
easy$error_sd^2     ## true
##plot: compare estimated eigenfunctions with the truth
eigenfest<-result$eigenfunctions
eigenf<-easy$eigenfunctions  ##true
par(mfrow=c(2,2))
for(i in 1:r){
  plot(grids.new,eigenfest[i,],ylim=range(eigenfest),xlab="time",ylab=paste("eigenfunction",i))
  points(grids, eigenf[i,],col=5,type="o")
}
##plot: compare estimated mean curve with the truth
muest<-result$fitted_mean
plot(grids.new,muest,xlab="time",ylab="mean curve",ylim=range(result$fitted_mean))
points(grids,numeric(length(grids)),col=5)
par(mfrow=c(1,1))
##look at the CV scores and convergence for each model: note that model (M=5, r=4) does not converge.
result$cv_scores   ##CV
result$converge   ##convergence
## derive fpc scores and look at the predicted curve
#fpc scores
fpcs<-fpca.score(easy$data,grids.new,muest,evalest,eigenfest,sig2est,r)
#get predicted trajectories on a fine grid: the same grid for which mean and eigenfunctions are evaluated
pred<-fpca.pred(fpcs, muest,eigenfest)
#get predicted trajectories on the observed measurement points
N<-length(grids.new)
par(mfrow=c(3,3))
for (i in 1:1){
  id<-i                                      ##for curve i
  t.c<-easy$data[easy$data[,1]==id,3]    ##measurement points
  t.proj<-ceiling(N*t.c)                     ##measurement points projected on the grid
  y.c<-easy$data[easy$data[,1]==id,2]    ##obs
  y.pred.proj<-pred[t.proj,id]               ##predicted obs on the measurement points
  #plots
  plot(t.c,y.c,xlim=c(0,1),ylim=range(pred[,id]),xlab="time",ylab="obs", main=paste("predicted trajectory of curve", id))
  points(grids.new,pred[,id],col=3,type='l'
  )
  ##points(t.c,y.pred.proj,col=2, pch=2)    ##predicted measurements at observed measurement times
}
par(mfrow=c(1,1))
############################## example II: "practical" case
##load data
## Not run:
data(prac)
##sample trajectory
plot(prac$data[,3],prac$data[,2],xlab="grid",ylab="",type='p')
for(i in 1:500){
  cur<-prac$data[prac$data[,1]==i,]
  points(cur[,3],cur[,2],type="l")
}
## candidate models for fitting
M.set<-c(5,10,15,20)
r.set<-5
##parameters for fpca.mle
ini.method="EM"
basis.method="bs"
sl.v=rep(0.5,10)
max.step=50
grid.l=seq(0,1,0.01)
grids=seq(0,1,0.002)
##fit candidate models by fpca.mle
result<-fpca.mle(prac$data, M.set,r.set,ini.method, basis.method,sl.v,max.step,grid.l,grids)
summary(result)
##rescaled grid
grids.new<-result$grid
##model selection result: the true model M=5, r=3 is selected with the smallest CV score among all converged models
M<-result$selected_model[1]
r<-result$selected_model[2]
##compare estimated eigenvalues with the truth
result$eigenvalues     ## estimated
prac$eigenvalues    ## true
##compare estimated error variance with the truth
result$error_var        ## estimated
prac$error_sd^2     ## true
##plot: compare estimated eigenfunctions with the truth
eigenf<-prac$eigenfunctions  ##true
par(mfrow=c(2,3))
for(i in 1:r){
  plot(grids.new,result$eigenfunctions[i,],ylim=range(result$eigenfunctions),xlab="time",ylab=paste("eigenfunction",i))
  points(grids, eigenf[i,],col=5,type="o")
}
##plot: compare estimated mean curve with the truth
plot(grids.new,result$fitted_mean,xlab="time",ylab="mean curve",ylim=range(result$fitted_mean))
points(grids,numeric(length(grids)),col=5)
par(mfrow=c(1,1))
##look at the CV scores and convergence for each model: note that model (M=5, r=4) does not converge.
result$cv_scores   ##CV
result$converge   ##convergence
## End(Not run)

#######################
library("fpca")
data("easy")
wide = fc.long2wide(easy$data[,1], easy$data[,3], easy$data[,2],bins = dgrid)
R = functionalImpute(wide, basis, K=3, maxIter = 50)

matplot(t(R$fit[10:15,]),t='l')
matplot(t(wide[10:15,]),t='p',pch = 'X',add=T)

plot(R$v[1,],t='l')
lines(R$v[2,])
lines(R$v[3,])
