sample_Us = function(Factors, genetic_effects, svd_ZZ_Ainv, Z_1){
  # %samples genetic effects (U) conditional on the factor scores F:
  #   % F_i = U_i + E_i, E_i~N(0,s2*(h2*ZAZ + (1-h2)*I)) for each latent trait i
  #   % U_i = zeros(r,1) if h2_i = 0
  #   % it is assumed that s2 = 1 because this scaling factor is absorbed in
  #   % Lambda
  #   % svd_ZZ_Ainv has parameters to diagonalize a*Z_1*Z_1' + b*I for fast
  #   % inversion:
  
  Q = svd_ZZ_Ainv$Q
  s1 = svd_ZZ_Ainv$s1
  s2 = svd_ZZ_Ainv$s2
  
  k = Factors$k
  n = genetic_effects$n
  tau_e = 1/(1-Factors$h2)
  tau_u = 1/Factors$h2
  b = t(Q) %*% Z_1 %*% t(Factors$scores * tau_e)
  z = matrix(rnorm(n*k), nrow=n)

  for (j in 1:k){
    if (tau_e[j] == 1){
      genetic_effects$U[j,] = matrix(rep(0,n), nrow=1)
    }
    else if (tau_e[j] == Inf){
      genetic_effects$U[j,] = Factors$scores[j,]
    }
    else {
      d = s2 %*% tau_u[j] + s1 %*% tau_e[j]
      mlam = b[,j] * (1/d)
      genetic_effects$U[j,] = t(Q %*% (mlam + z[,j] * (1/sqrt(d))))
    }
  }
  
  return(genetic_effects)
}