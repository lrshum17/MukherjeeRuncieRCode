sample_means = function(Ytil, Qt_Design, N, resid, random_precision, svd_Design_Ainv){
  # %when used to sample [B;D]:
  #   %Y - FL' - Z_2W = XB + ZD + E, vec(E)~N(0,kron(Psi_E,In)). 
  # % Note: conditioning on F, L and W.
  # %The vector [b_j;d_j] is sampled simultaneously. Each trait is sampled separately because their
  # %conditional posteriors factor into independent MVNs.
  # %note:svd_Design_Ainv has parameters to diagonalize mixed model equations for fast inversion: 
  # %inv(a*blkdiag(fixed_effects.cov,Ainv) + b*[X; Z_1][X; Z_1]') = Q*diag(1./(a.*s1+b.*s2))*Q'
  # %Qt_Design = Q'*Design, which doesn't change each iteration. Design = [X;Z_1]
  # %
  # %function also used to sample W:
  # %Y - FL' - XB - ZD = Z_2W + E, vec(E)~N(0,kron(Psi_E,In)). 
  # %Here, conditioning is on B and D.
  p=resid$p
  Q = svd_Design_Ainv$Q
  s1 = svd_Design_Ainv$s1
  s2 = svd_Design_Ainv$s2
  
  means = (Qt_Design %*% t(Ytil)) * t(resid.ps)
  location_sample = matrix(rep(0,N*p), nrow = N)
  Zlams = matrix(rnorm(N*p), nrow = N)
  for (j in 1:p){
    d = s1 %*% random_precision[j] + s2 %*% resid$ps[j]
    mlam = means[,j] * (1/d)
    location_sample[,j] = Q %*% (mlam + Zlams[,j] * (1/sqrt(d)))
  }
  location_sample = t(location_sample)
  return(location_sample)
}