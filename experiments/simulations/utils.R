dyadic.basis = function(A, YA, min.n = 200) {
  indicator = matrix(1, nrow(A), 1)

  # stop splitting
  if (nrow(A) <= min.n / 2) {
    return(indicator)
  } else if (nrow(A) <= min.n) {
    # just one more greedy split
    improvement = sapply(1:ncol(A), function(j) {
      med.j = median(A[,j])
      if(max(A[,j]) == med.j) return (0)
      (mean(YA[A[,j] <= med.j]) - mean(YA[A[,j] > med.j]))^2
    })
    if (max(improvement) == 0) {
      return(indicator)
    }
    # do the greedy split
    j.max = which.max(improvement)
    med.j.max = median(A[,j.max])
    is.left = (A[,j.max] <= med.j.max)
    return(cbind(indicator, as.numeric(is.left), as.numeric(!is.left)))
  }
  sub.basis = lapply(1:ncol(A), function(j) {
    med.j = median(A[,j])
    if(max(A[,j]) == med.j) return (rep(1, nrow(A)))
    is.left = (A[,j] <= med.j)
    left = dyadic.basis(A[is.left,], YA[is.left], min.n)
    right = dyadic.basis(A[!is.left,], YA[is.left],  min.n)
    all = matrix(0, nrow(A), ncol(left) + ncol(right))
    all[is.left, 1:ncol(left)] = left
    all[!is.left, ncol(left) + 1:ncol(right)] = right
    all
  })
  cbind(indicator, Reduce(cbind, sub.basis))
}

make.tree.basis = function(X, W, num.tree = 5) {
  nobs = length(W)

  Reduce(cbind, lapply(1:num.tree, function(rep) {
    if (rep == 1) {
      features = 1:ncol(X)
    } else {
      features=sample(1:ncol(X), ceiling(3 * ncol(X)/4))
    }

    DF=data.frame(W=W, X[,features])
    tree = rpart(W~., data = DF, control=rpart.control(cp=0))
    nodes = matrix(0, nobs, nrow(tree$frame))
    for(iter in 1:nobs) {
      nodes[iter, tree$where[iter]] = 1
    }

    node.names = as.numeric(rownames(tree$frame))
    for(idx in nrow(tree$frame):2) {
      parent.name = floor(node.names[idx]/2)
      parent.idx = which(node.names == parent.name)
      nodes[,parent.idx] = nodes[,parent.idx] + nodes[,idx]
    }
    nodes
  }))
}

# helper function to run simulations
diffindiff.comparison = function(X, Y, Ti, Si, f, h, b, TAU, tau,
                          P_T1_S1, P_T1_S0, P_T0_S1, P_T0_S0, alg,
                          order=NULL) {
  if(is.null(order)) {
    order = 1
    dd = min(1000, nrow(X) - ncol(X))
    while (dd > 0 & order < 10) {
      order = order + 1
      dd = dd - choose(ncol(X) - 1 + order, order)
    }
  }
  # create oracle weights
  gamma.oracle = (Ti*Si)/P_T1_S1 + (Ti*(1-Si))/P_T1_S0 + ((1-Ti)*Si)/P_T0_S1 + ((1-Ti)*(1-Si))/P_T0_S0
  Basis = generate.basis(X,order)
  Basis = cbind(1,Basis)

  P_T1 = P_T1_S0 + P_T1_S1
  P_T0 = P_T0_S0 + P_T0_S1
  P_S1 = P_T1_S1 + P_T0_S1
  P_S0 = P_T1_S0 + P_T0_S0
  if (alg == "oracle_amle"){
    oracle_amle = oracle_learner(X=Basis, Y=Y, Ti=Ti, Si=Si, f=f, h=h, b=b, TAU=TAU, tau=tau, gamma=gamma.oracle,
                                 P_T1_S1=P_T1_S1, P_T1_S0=P_T1_S0, P_T0_S1=P_T0_S1, P_T0_S0=P_T0_S0,
                                 constant_eff = "non_constant")
    return(oracle_amle)
  } else if (alg == "IPW"){
    gamma.plugin = plugin(Basis, Ti, Si)
    return(list(TAU_hat = mean(gamma.plugin*Y), std.err.est = NA))
  }  else if (alg == "DiD_AIPW") {
    gamma.plugin = plugin(Basis, Ti, Si)
    DiD_AIPW = DiD(X=Basis, Y=Y, Ti=Ti, Si=Si,
                                     constant_eff="non_constant", gamma=gamma.plugin)
    return(DiD_AIPW)
  }  else if (alg == "DiD_AMLE") {
    gamma.minimax = minimax(Basis, Ti, Si)
    DiD_AMLE = DiD(X=Basis, Y=Y, Ti=Ti, Si=Si,
                                     constant_eff="non_constant", gamma=gamma.minimax)
    return(DiD_AMLE)
  } else {
    stop("bad alg for non-const")
  }
}
