require(pracma)

save_posterior_samples = function(sp_num, params, Posterior, resid, fixed_effects,
                                                       genetic_effects, Factors, interaction_effects){
  #%save posteriors. Re-scale samples back to original variances.
  
  #%save factors
  sp = params$sp
  VY = params$VY       
  Lambda = Factors$Lambda * sqrt(t(VY))     #%re-scale by Y variances
  G_h2 = Factors$h2
  U = genetic_effects$U     
  delta = Factors$delta
  genetic_ps = genetic_effects$ps/t(VY)
  resid_ps = resid$ps/t(VY)
  
  #%save factor samples
  Lambda = Lambda[, 1:Factors$k]
  if(prod(dim(Lambda)) > nrow(Posterior$Lambda)){
    #%expand factor sample matrix if necessary
    Posterior$Lambda = rbind(Posterior$Lambda, zeros(prod(dim(Lambda))-nrow(Posterior$Lambda), sp))
    Posterior$U = rbind(Posterior$U, zeros(prod(dim(U))-nrow(Posterior$U), sp))
    Posterior$delta = rbind(Posterior$delta, zeros(prod(dim(delta))-nrow(Posterior$delta), sp))
    Posterior$G_h2 = rbind(Posterior$G_h2, zeros(prod(dim(G_h2))-nrow(Posterior$G_h2), sp))
  }
  
  Posterior$Lambda[1:prod(dim(Lambda)), sp_num] = Lambda(:) # CHECK!!!!!!
  Posterior$delta[1:prod(dim(delta)), sp_num] = delta
  Posterior$G_h2[1:prod(dim(G_h2)), sp_num] = G_h2
  Posterior$U[1:prod(dim(U)), sp_num] = U(:) # check!!!!!
  
  Posterior$ps[,sp_num] = genetic_ps
  Posterior$resid_ps[,sp_num] = resid_ps
  
  #%save B,U,W
  Posterior$B = (Posterior$B * (sp_num-1) + (fixed_effects$B * sqrt(t(VY))))/sp_num
  Posterior$d = (Posterior$d * (sp_num-1) + (genetic_effects$d * sqrt(t(VY))))/sp_num
  Posterior$W = (Posterior$W * (sp_num-1) + (interaction_effects$W * sqrt(t(VY))))/sp_num
  
  if(sp_num %% 100 == 0){
    save('Posterior', 'Posterior', 'params') #CHECK
  }

  return(Posterior)
}