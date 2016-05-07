sample_factors_scores = function(Ytil, Factors, resid, genetic_effects, Z_1){
  #Sample factor scores given factor loadings, U, factor heritabilities and
  #phenotype residuals
  k = Factors$k
  n = Factors$n
  Lambda = Factors$Lambda
  Lmsg = Lambda * resid$ps
  tau_e = 1/(1-Factors$h2)
  S = t(chol(t(Lambda) %*% Lmsg+diag(tau_e)))
  Meta = solve(t(S)) %*% (solve(S) %*% ((t(Lmsg) %*% Ytil + genetic_effects$U * Z_1 * tau_e)))
  Factors$scores = Meta + solve(t(S)) %*% matrix(rnorm(k*n), nrow = k) 
  return(Factors)
}