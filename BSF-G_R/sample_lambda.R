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
    vlam = inv(Llam) %*% means
    mlam = inv(t(Llam)) %*% vlam
    ylam = inv(t(Llam)) %*% Zlams[,j]
    Factors$Lambda[j,] = (ylam + mlam)
  }
}










# 
# function [Factors] = sample_lambda( Ytil,Factors, resid,genetic_effects,eig_ZAZ )
# %Sample factor loadings (Factors.Lambda) while marginalizing over residual
# %genetic effects: Y - Z_2W = FL' + E, vec(E)~N(0,kron(Psi_E,In) + kron(Psi_U, ZAZ^T))
# %note: conditioning on W, but marginalizing over U.
# %sampling is done separately by trait because each column of Lambda is
# %independent in the conditional posterior
# %note: eig_ZAZ has parameters that diagonalize aI + bZAZ for fast
# %inversion: inv(aI + bZAZ) = 1/b*Ur*diag(1./(eta+a/b))*Ur'
# 
# p=resid.p;
# k=Factors.k;
# 
# Ur = eig_ZAZ.vectors;
# eta = eig_ZAZ.values;
# FtU = Factors.scores*Ur;
# UtY = Ur'*Ytil';
# 
# Zlams = normrnd(0,1,k,p);
# for j = 1:p,
# FUDi = genetic_effects.ps(j)*bsxfun(@times,FtU,1./(eta' + genetic_effects.ps(j)/resid.ps(j)));
#                                                    means = FUDi*UtY(:,j);
#                                                    Qlam = FUDi*FtU' + diag(Factors.Plam(j,:)); 
#                                                    Llam = chol(Qlam,'lower');
#                                                    vlam = Llam\means; mlam = Llam'\vlam; ylam = Llam'\Zlams(:,j);
#                                                    Factors.Lambda(j,:) = (ylam + mlam);
#                                                    end
#                                                    end
#                                                    