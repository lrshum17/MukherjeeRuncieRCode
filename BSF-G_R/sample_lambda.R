sample_lambda = function(Ytil, Factors, resid, genetic_effects, eig_ZAZ){
  p=resid$p
  k=Factors$k
  
  Ur = eig_ZAZ$vectors
  eta = eig_ZAZ$values
  FtU = Factors$scores*Ur
  UtY = t(Ur) %*% t(Ytil)
  
  Zlams = matrix(rnorm(k*p), nrow = k);
  for (j in 1:p){
    FUDi = genetic_effects$ps[j] %*% (FtU * (1/(t(eta) + genetic_effects$ps[j]/resid$ps[j])))
    means = FUDi %*% UtY[,j]
    Qlam = FUDi %*% t(FtU) + diag(Factors$Plam[j,])
    Llam = t(chol(Qlam))
    vlam = solve(Llam) %*% means
    mlam = solve(t(Llam)) %*% vlam
    ylam = solve(t(Llam)) %*% Zlams[,j]
    Factors$Lambda[j,] = (ylam + mlam)
  }
  
  return(Factors)
}