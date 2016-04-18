#Need to define Z_1 before using

sample_factors_scores = function(Ytil, Factors, resid, genetic_effects, Z_1){
  #Sample factor scores given factor loadings, U, factor heritabilities and
  #phenotype residuals
  k = Factors$k
  n = Factors$n
  Lambda = Factors$Lambda
  Lmsg = Lambda * resid$ps
  tau_e = 1/(1-Factors$h2)
  S = t(chol(t(Lambda) %*% Lmsg+diag(tau_e)))
  Meta = inv(t(S)) %*% (inv(S) %*% ((t(Lmsg) %*% Ytil + genetic_effects$U * Z_1 * tau_e)))
  Factors$scores = Meta + inv(t(S)) %*% matrix(rnorm(k*n), nrow = k) 
  return(Factors)
}

# function [ Factors ] = sample_factors_scores( Ytil, Factors,resid,genetic_effects )
# %Sample factor scores given factor loadings, U, factor heritabilities and
# %phenotype residuals
# global Z_1
# k = Factors.k;
# n = Factors.n;
# Lambda = Factors.Lambda;
# Lmsg = bsxfun(@times,Lambda,resid.ps);
# tau_e = 1./(1-Factors.h2);
# S=chol(Lambda'*Lmsg+diag(tau_e),'lower');
# Meta = S'\(S\(Lmsg'*Ytil + bsxfun(@times,genetic_effects.U*Z_1,tau_e)));
# Factors.scores = Meta + S'\randn(k,n);   
# end
#               