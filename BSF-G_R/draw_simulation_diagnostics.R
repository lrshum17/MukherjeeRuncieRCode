draw_simulation_diagnostics = function(sp_num,params,Factors,genetic_effects,resid,
                                       Posterior,gen_factor_Lambda,error_factor_Lambda,G,R,h2){
  #draw some diagnostic plots   
  VY = params$VY
  p = params$p
  Lambda = Factors$Lambda
  G_h2s = t(Factors$h2)
  E_h2s = 1-G_h2s
  G_Lambda = Factors$Lambda * sqrt(G_h2s)
  E_Lambda = Factors$Lambda * sqrt(E_h2s)
  actual_G_Lambda = gen_factor_Lambda * 1/sqrt(t(VY))
  actual_E_Lambda = error_factor_Lambda * 1/sqrt(t(VY))
  
  G_est = G_Lambda %*% t(G_Lambda) + diag(1/genetic_effects$ps)
  G_act = G
  G_act = G_act/sqrt(t(VY) %*% VY)
  
  E_est = E_Lambda %*% t(E_Lambda) + diag(1/resid$ps)
  E_act = R
  E_act = E_act/sqrt(t(VY) %*% VY)
  
  delta=Factors$delta
  
  #Figure 1: accuracy of current parameter estimates
  #plot 1: correlation between each estimated factor and the true factors.
  #Outer row and column shows maximum correlation for each fitted or actual
  #factor, respectively.
  jpeg('Diagnostics_Lambda.jpg')
  p_view = 1:p
  f1_row=3
  f1_col=3
  par(mfrow = c(f1_row,f1_col))
  
  clims=c(0, 1)
  cors = abs(Cor(t(Lambda),t(actual_E_Lambda)))
  cors = rbind(cors,t(apply(cors,2,max)))
  cors = cbind(cors,t(t(apply(cors,1,max))))
  image(cors,zlim = clims)
  #Check colorbar
  
  #plot 2: visualize estimated genetic covariance matrix (as correlation matrix)
  clims=c(-1, 1)
  image(CovToCor(G_est),zlim = clims)
  
  #plot 3: visualize estimated residual covariance matrix (as correlation matrix)
  clims=c(-1, 1)
  image(CovToCor(E_est),zlim = clims)
  
  #plot 4: visualize estimated genetic factor loading matrix
  clim=max(0.1,min(.75,max(abs(G_Lambda))));
  clims=c(-clim, clim)
  image(G_Lambda[p_view,],zlim = clims)
  #colorbar Check
  
  #plot 5: visualize actual genetic factor loading matrix
  image(cbind(actual_G_Lambda[p_view,], matrix(0, nrow=length(p_view), ncol=1)),zlim = clims)
  #colorbar Check
  
  #plot 6: visualize estimated residual factor loading matrix
  clim=max(0.1,min(.75,max(abs(E_Lambda))))
  clims=c(-clim, clim)
  image(E_Lambda[p_view,],zlim = clims)
  #colorbar Check
  
  #plot 7: Element-wise comparison of estimated and actual G matrices
  matplot(G_est,G_act,pch=20)
  abline(a=0, b=1)
  
  #plot 8: Element-wise comparison of estiamted and actual R matrices
  matplot((E_ect+G_ect)*upper.tri(E_ect+G_ect),(E_act+G_act)*upper.tri(E_act+G_act),pch=20)
  abline(a=0, b=1)
  
  #plot 9: Plot of factor heritabilities
  matplot(G_h2s, ty="l", col='red', ylim = c(0, 1))
  matlines(E_h2s,col='blue')
  dev.off()
  
  if (sp_num <= 1){
    k = min(ncol(actual_E_Lambda),ncol(Lambda))
    max_cors = matrix(0,nrow=k,ncol=1)
    cors = cors[1:nrow(cors)-1,1:ncol(cors)-1]
    for (j in 1:k){
      max_cors[j] = which.max(cors[j,])
    }
    o = order(max_cors)
    # figure(5) # Check
    # plot(Lambda[,o], actual_E_Lambda,pch=20)
    L = Lambda[,o]
    
    #figure(6) Check
    for (j in 1:min(k,5*5)){
      subplot(5,5,j) # Check
      plot(actual_E_Lambda[,j],L[,j],pch=20)
      abline(a=0, b=1)
      abline(a=0, b=-1)
      end
    }
  } 
  
  if (sp_num>1){
    #Figure of trace plots of the largest elements of each column of the factor loading matrix
    jpeg('Diagnostics_autoCorr_Lambda.jpg', width = 600, height = 800)
    f2_row=8
    f2_col=2
    par(mfrow=c(f2_row,f2_col))
    for (k in 0:(min(2*f2_row,nrow(Posterior$Lambda)/p)-1)){
      o = order(-abs(apply(Posterior$Lambda[(1:p)+k*p,max(1,sp_num-100):sp_num],1,mean)))
      matplot(t(Posterior$Lambda[o[1:5]+k*p,1:sp_num]),ty="l")
    }
    dev.off()
    
    #Figure of trace plots and histograms of the factor heritablities
    jpeg('Diagnostics_factor_selection.jpg', width = 600, height = 800)
    f4_row=8
    f4_col=4
    par(mfrow=c(f4_row,f4_col))
    for (k in 1:min(2*f4_row,ncol(Posterior$G_h2))){
      h2s = Posterior$G_h2[k,1:sp_num]
      if (sum(h2s[!is.na(h2s)])==0){
        next()
      }
      plot(h2s, ty="l", ylim=c(0, 1))
      hist(h2s, breaks = 100, xlim=c(-0.1, 1))
    }
    dev.off()
    
    #Figure similar to figure 1, but posterior mean parameter estimates 
    jpeg('Diagnostics_posterior_mean.jpg')
    f1_row=3;
    f1_col=2;
    
    k = ncol(Posterior$Lambda)/p
    h2s = Posterior$G_h2[,1:sp_num] 
    G_Lambdas = matrix(0,nrow = nrow(Posterior$Lambda), ncol = ncol(Posterior$Lambda))
    Lambda_est = matrix(0,p,k)
    G_est = matrix(0,p,p)
    E_est = matrix(0,p,p)
    
    for (j in 1:sp_num){
      Lj = matrix(Posterior$Lambda[,j],p,k) * (1/sqrt(t(VY)))
      h2j = Posterior$G_h2[,j]
      G_Lj = Lj %*% diag(sqrt(h2j))
      G_Lambdas[,j] = matrix(G_Lj,nrow=nrow(G_Lj)*ncol(G_Lj))
      Gj = G_Lj %*% t(G_Lj) + diag(1/(t(VY) * Posterior$ps[,j]))
      G_est = G_est + Gj/sp_num
      
      E_Lj = Lj %*% diag(1-h2j) %*% t(Lj) + diag(1/(t(VY) * Posterior$resid_ps[,j]))
      E_est = E_est + E_Lj/sp_num
      Lambda_est = Lambda_est + matrix(Posterior$Lambda[,j],p,k)/sp_num
    }
    # G_est = 4 * G_est
    G_Lambda = matrix(apply(G_Lambdas,1,mean),p,k)
    # Lambda = Lambda * sqrt(t(VY))
    
    #plot 1: visualize estimated genetic factor loading matrix
    clims=c(-1, 1)
    image(CovToCor(G_est),zlim = clims)
    
    #plot 2: visualize actual genetic factor loading matrix
    image(cbind(actual_G_Lambda[p_view,], rep(0,length(p_view))), zlim = clims)
    # Check colorbar;
    
    #plot 3: Element-wise comparison of estiamted and actual G matrices
    matplot((G_ect)*upper.tri(G_ect),(G_act)*upper.tri(G_act),pch=20)
    abline(a=0, b=1)  
    
    #plot 4: visualize estimated genetic factor loading matrix
    clim=min(.75,max(abs(G_Lambda)))
    clims=c(-clim, clim)
    image(G_Lambda[p_view,], zlim = clims)
    # Check lx = xlim;
    # Check colorbar;  

    #plot 5: Plot of trait heritabilities
    plot(diag(G_est)/(diag(G_est)+diag(E_est)),h2, pch=20,col='red')
    abline(a=0, b=1) 
    
    #plot 6: Plot of factor heritabilities
    plot(apply(h2s,1,mean),ty="l",col = "red", ylim = c(0, 1))
    lines(1-(apply(h2s,1,mean)),col = "blue")
    # Check colorbar
    # Check xlim(lx)
    
    dev.off()
    
    clims=c(0, 1)        
    cors = abs(Cor(t(Lambda_est),t(error_factor_Lambda)))
    cors = rbind(cors, apply(cors,2,max))
    cors = cbind(cors, t(apply(t(cors),2,max)))
    
    k=min(ncol(actual_E_Lambda),ncol(Lambda_est))
    max_cors = matrix(0,k,1)
    cors = cors[1:(nrow(cors)-1),1:(ncol(cors)-1)]
    for (j in 1:k){
      max_cors[j] = which.max(cors[,j])
    }
    o = order(max_cors)
    # figure(5)
    # plot(Lambda[,o],actual_E_Lambda,pch=20)
    L=error_factor_Lambda[,o]
    # a= solve(t(actual_E_Lambda) %*% actual_E_Lambda) %*% t(actual_E_Lambda) %*% Lambda
    # plot(apply(abs(a),2,max))
    
    #figure(6)
    for (j in 1:min(k,5*5)){
      plot(L[,j],Lambda_est[,j],pch=20)
      abline(a=0,b=1)
      abline(a=0,b=-1)
    }
  }
}