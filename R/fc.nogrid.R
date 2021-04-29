library("fda")
library("tidyr")

gradient.clean <- function(B,W,y,index){
  ## B is M x K matrix of basis evaluations, at the M = \sum_{i=1}^N n_i timepoints.
  ## W is N x K current solution matrix
  ## y is M vector of responses (measurements)
  ## index is M vector, indicating which elements belong to which "subject" ie 1,1,1,2,2,3,3,3 ..., N,N,N
  resid <- y - drop( rowSums(B * W[index, ]))
  rB <- resid * B # elementwise
  lrB <- split(rB, index) # drops the dimensions
  K <- ncol(B)
  G <- sapply(lrB, function(x)colSums(matrix(x,ncol=K)))
  list(grad=t(G), res=sum(resid**2))
}

nogrid.fimpute.fit = function(data,
                       value.vars,
                       time.var,
                       id.var,
                       basis = NULL,
                       niter = 50,
                       pp = 10,
                       lr = 0.05,
                       lambda = 1,
                       tol = 1e-7,
                       dgrid = 100,
                       nogrid = TRUE,
                       verbose = FALSE
){
#  cat(value.vars, time.var, id.var, tol)

  # Prepare a grid and a basis
  rangeval = c(min(data[[time.var]],na.rm=TRUE),max(data[[time.var]],na.rm=TRUE))
  grid = seq(rangeval[1], rangeval[2], length.out = dgrid)

  if (is.null(basis)){
    basis = create.bspline.basis(rangeval = rangeval, nbasis=7)
  }

  basis.on.grid = eval.basis(grid, basis)
  df = basis$nbasis

  # Estimate mean curves and subtract
  mean_curves = list()
  for (value.var in value.vars){
    nna = !is.na(data[[value.var]])

    # Spline version
    mean_curves[[value.var]] = smooth.spline(data[[time.var]][nna], data[[value.var]][nna], df=2)
    data[[value.var]] = data[[value.var]] - predict(mean_curves[[value.var]], data[[time.var]])$y
  }

  # Convert to long, sort and reindex
  data_long = gather(data, measurement, value, value.vars, factor_key=TRUE)
  data_long = na.omit(data_long)
  data_long = data_long[order(data_long[[id.var]], data_long[[time.var]]),]
  data_long$idx = as.numeric(factor(data_long[[id.var]]))

  # Count remaining subjects
  subjects = unique(data_long$idx)
  n = length(subjects)
  nvars = length(value.vars)
  pp = min(pp, df*nvars)

  # Store solution matrix
  Wold = matrix(rnorm(n*df*nvars)/100,n,df*nvars)

  # Evaluate basis on the observe timpoints and normalize variables
  basis_evals = list()
  slices = list()
  mm = list()
  ssd = list()
  Y = list()

  for (value.var in value.vars){
    slices[[value.var]] = data_long[["measurement"]] == value.var

    mm[[value.var]] = mean(data_long[["value"]][slices[[value.var]]])
    ssd[[value.var]] = sd(data_long[["value"]][slices[[value.var]]])
    data_long[["value"]][slices[[value.var]]] = (data_long[["value"]][slices[[value.var]]] - mm[[value.var]])/ssd[[value.var]]

    #    basis_evals[[value.var]] = predict(basis,
    #                                       data_long[[time.var]][slices[[value.var]]])

    if(nogrid == TRUE){
      basis_evals[[value.var]] = eval.basis(data_long[[time.var]][slices[[value.var]]], basis)
    }
    else {
      basis_evals[[value.var]] = eval.basis(grid, basis)
#      M = data.frame(idx = data_long$idx[slices[[value.var]]],
#        vals = data_long$value[slices[[value.var]]],
#        time = sapply(data_long$time[slices[[value.var]]], function(x){ which.min(abs(x-grid)) }))
      Y[[value.var]] = fcomplete:::fc.long2wide(data_long$idx[slices[[value.var]]], data_long$time[slices[[value.var]]], data_long$value[slices[[value.var]]], bins = length(grid))
    }
  }

  ## Estimate W directly (no grid)
  for (i in 1:niter){
    total_grad = Wold
    total_grad[] = 0
    total_res = 0

    for (j in 1:nvars){
      value.var = value.vars[j]
      if (nogrid == TRUE){
        res.update = gradient.clean(basis_evals[[value.var]],
                                   Wold[,((j-1)*df+1):(j*df)],
                                   data_long[["value"]][slices[[value.var]]],
                                   data_long$idx[slices[[value.var]]])
      }
      else {
        # Get current estimate using W
        Yhat = Wold[,((j-1)*df+1):(j*df)] %*% t(basis_evals[[value.var]])
        yobs = !is.na(Y[[value.var]])
        Yres = Yhat
        Yres[] = 0
        Yres[yobs] = Y[[value.var]][yobs] - Yhat[yobs]
        res.update = list()
        res.update$grad = fcomplete:::project.on.basis(Yres, basis_evals[[value.var]])
        res.update$res = sum(Yres) / length(yobs)
      }
      total_grad[unique(data_long$idx[slices[[value.var]]]),((j-1)*df+1):(j*df)] = res.update$grad
      total_res = total_res + res.update$res
    }
    W = Wold + total_grad*lr

    if (i %% 1 == 0 && verbose){
      print(paste("Iter",i,"Fit:", total_res, ", obj:", total_res + lambda*norm(W)))
    }

    # Truncated SVD
    ss = svd(W, nu = pp, nv = pp)
    ss$d = ss$d - lambda*lr
    ss$d[ss$d<0] = 0

    dd = diag(ss$d[1:pp])
    if (length(ss$d[1:pp]) == 1)
      dd = ss$d[1:pp]

    W = as.matrix(ss$u) %*% dd %*% t(as.matrix(ss$v))

    if ( sum((W-Wold)**2) <= tol * sum(Wold**2))
      break
    Wold = W
  }


  model_functions = list()
  for (j in 1:nvars){
    value.var = value.vars[j]
    model_functions[[value.var]] = t(t(Wold[,((j-1)*df+1):(j*df)] %*% t(basis.on.grid) * ssd[[value.var]] + mm[[value.var]]) + predict(mean_curves[[value.var]], x=grid)$y)
    data_long[["value"]][slices[[value.var]]] = data_long[["value"]][slices[[value.var]]]*ssd[[value.var]] + mm[[value.var]]+
      predict(mean_curves[[value.var]], data_long[[time.var]][slices[[value.var]]])$y
  }
  # for compatibility
  V = svd(model_functions[[1]])$v
  list(u=Wold,
       v=V,
       fit=model_functions,
       mean_curves=mean_curves,
       data_long=data_long,
       subjects=subjects,
       slices=slices,
       basis=basis,
       grid=grid
  )
}

# TODO: optimize
predict.fimpute = function(fit,grid,newdata,time.var,id.var){
  preds = c()
  for (i in 1:nrow(newdata)){
    g = which(grid > newdata[i,time.var])[1]
    if (is.null(g) || is.na(g)){
      g = 1
    }
    else if(g > ncol(fit)){
      g = ncol(fit)
    }
    preds = c(preds, fit[newdata[i,id.var],g])
  }
  preds
}

cv.nogrid.fimpute = function(data,
                      value.vars,
                      time.var,
                      id.var,
                      ...,
                      lambdas=c(0,1),
                      val.ratio=0.05,
                      nreps = 5){
  bestModel = NULL
  bestL = 1e10
  loss = c()
  #  functions.on.grid = t(eval.fd(0:99/99,functions))

  nsubj = length(unique(data[[id.var]]))
  test.masks = list()

  for (i in 1:nreps){
    test.subjects = 1:nsubj %in% sample(1:nsubj)[1:floor(nsubj*val.ratio)]
    test.idx = c()
    for (subj in test.subjects){
      test.idx = c(test.idx, sample(which(data[[id.var]] == subj))[1])
    }
    test.masks[[i]] = 1:length(data[[id.var]]) %in% test.idx
  }

  for (lambda in lambdas){
    l = 0
    for (i in 1:nreps){
      test.mask = test.masks[[i]]
      model = nogrid.fimpute.fit(data[!test.mask,],
                          value.vars,
                          time.var,
                          id.var,
                          ...,
                          lambda=lambda)

      # TODO: for each var
      for (value.var in value.vars){
        l = l + sum((data[test.mask,] - predict.fimpute(model$fit[[value.var]], model$grid, data[test.mask,], time.var, id.var))**2,na.rm=TRUE)
      }
    }
    l = l / nreps

    ## For tests when ground truth is known
    # l = mean((functions.on.grid - model$functions$measurement)**2) / mean((functions.on.grid - 0)**2)
    loss = c(loss, l)
    if (l < bestL){
      bestL = l
      bestLambda = lambda
#      bestModel = model
    }
  }
  model = nogrid.fimpute.fit(data,
                             value.vars,
                             time.var,
                             id.var,
                             ...,
                             lambda=bestLambda)

  list(model=model, loss=loss)
}

