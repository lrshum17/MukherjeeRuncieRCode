library(graphics)

draw_results_diagnostics = function(sp_num,params,Factors,Posterior){
  VY = params$VY
  p = params$p
  f1_row=3
  f1_col=2
  p_view = 1:p
  
  Lambda = Factors$Lambda
  G_h2s = t(Factors$h2)
  E_h2s = 1-t(Factors$h2)
  G_Lambda = Factors$Lambda*sqrt(G_h2s)
  E_Lambda = Factors$Lambda*sqrt(E_h2s)
  delta = Factors$delta
  
  
  
  if(sp_num > 10){
    #Figure of trace plots of the largest elements of each column of the factor loading matrix
    f2_row = 8
    f2_col = 2
    jpeg('Diagnostics_autoCorr_Lambda.jpg', width = 600, height = 800)
    par(mfrow=c(8,2))
    for(k in 0:(min(2*f2_row,nrow(Posterior$Lambda)/p)-1)){
      o = order(-abs(apply(Posterior$Lambda[(1:p)+k*p,max(1,sp_num-100):sp_num],2,mean)))
      matplot(Posterior$Lambda[o[1:5]+k*p,1:sp_num], ty="l")
    }
    dev.off()
    
    
    
    #Figure of trace plots and histograms of the factor heritablities
    jpeg('Diagnostics_factor_selection.jpg', width = 600, height = 800)
    par(mfrow=c(8,4))
    for (k in 1:min(2*f4_row,nrow(Posterior$G_h2))){
      h2s = Posterior$G_h2[k,1:sp_num]
      if (sum(h2s(!is.na(h2s)))==0){
        next
      }
      plot(h2s, ty="l", ylim = c(0,1))
      hist(h2s, breaks = 100, xlim = c(0,1))
    }
    dev.off()
    
    
    k = nrow(Posterior$Lambda)/p
    h2s = Posterior$G_h2[,1:sp_num]
    G_Lambdas = matrix(0, nrow = nrow(Posterior$Lambda), ncol = ncol(Posterior$Lambda))
    G_est = matrix(rep(0,p*p), nrow = p)
    E_est = matrix(rep(0,p*p), nrow = p)
    for (j in 1:sp_num){
      Lj = matrix(Posterior$Lambda[,j],nrow = p, ncol = k)*(1/sqrt(t(VY)))
      h2j = Posterior$G_h2[,j]
      G_Lj = Lj %*% diag(sqrt(h2j))
      G_Lambdas[,j] = matrix(G_Lj, nrow = nrow(G_Lj)*ncol(G_Lj))
      Gj = G_Lj %*% t(G_Lj) + diag(1/(t(VY) * Posterior$ps[,j]))
      G_est = G_est + Gj/sp_num
      
      E_Lj = Lj %*% diag(1-h2j) %*% t(Lj) + diag(1/(t(VY) * Posterior$resid_ps[,j]))
      E_est = E_est + E_Lj/sp_num
    }
    G_Lambda = matrix(apply(G_Lambdas, 2, mean), nrow = p, ncol = k)
    
  }
  
  
  
  
  #Figure 1: accuracy of current parameter estimates
  #plot 1: estimated number of important factors
  jpeg('Diagnostics_Lambda.jpg')
  par(mfrow = c(f1_row,f1_col))
  plot(Factors$no_f, ty="l")
  
  if(sp_num > 10){
    #plot 2: visualize posterior mean genetic correlation matrix
    clims=c(-1, 1)
    image(CovToCor(G_est), zlim = clims)
    #Need to laod CovToCor.R
  }
  
  
  #plot 3: visualize estimated genetic factor loading matrix
  clim = max(.1,min(.75,max(max(abs(G_Lambda)))))
  clims = c(-clim, clim)
  image(G_Lambda[p_view,], zlim = clims)
  #Lack of a colorbar
  
  if(sp_num > 10){
    #plot 4: visualize posterior mean genetic factor loading matrix
    clim = min(.75,max(max(abs(G_Lambda))))
    clims = c(-clim, clim)
    image(G_Lambda[p_view,], zlim = clims)
    #lx = xlim;
    #Lack of a colorbar
  }
  
  #plot 5: visualize estimated residual factor loading matrix
  clim = max(.1,min(.75,max(max(abs(E_Lambda)))))
  clims = c(-clim, clim)
  image(E_Lambda[p_view,], zlim = clims)
  #Lack of a colorbar
  
  #plot 6: Plot of factor heritabilities
  plot(G_h2s, ty="l", col="red", ylim = c(0, 1))
  lines(E_h2s, col='blue')
  
  if(sp_num > 10){
    #plot 6: posterior mean factor heritabilities
    plot(apply(h2s,2,mean),ty="l",col = 'red', ylim=c(0, 1))
    lines(1-apply(h2s,2,mean),col = 'blue')
    #colorbar
    #xlim(lx)
  }
  dev.off()

  
}       
  
  
  
  
  









# 
# 
# function draw_results_diagnostics(~,sp_num,params,Factors,Posterior)
# 
# VY = params.VY;
# p = params.p;
# fig1=figure(1);
# f1_row=3;
# f1_col=2;
# clf
# p_view = 1:p;
# 
# Lambda = Factors.Lambda;
# G_h2s = Factors.h2';
# E_h2s = 1-Factors.h2';
# G_Lambda = bsxfun(@times,Factors.Lambda,sqrt(G_h2s));
# E_Lambda = bsxfun(@times,Factors.Lambda,sqrt(E_h2s));
# delta=Factors.delta;
# 
# %Figure 1: accuracy of current parameter estimates
# %plot 1: estimated number of important factors
# subplot(f1_row,f1_col,1)
# plot(Factors.no_f)
# 
# %plot 3: visualize estimated genetic factor loading matrix
# subplot(f1_row,f1_col,3)
# clim=max(.1,min(.75,max(max(abs(G_Lambda)))));
# clims=[-clim clim];
# imagesc(G_Lambda(p_view,:),clims)
# colorbar;
# 
# %plot 5: visualize estimated residual factor loading matrix
# subplot(f1_row,f1_col,5)
# clim=max(.1,min(.75,max(max(abs(E_Lambda)))));
# clims=[-clim clim];
# imagesc(E_Lambda(p_view,:),clims)
# colorbar;
# 
# %plot 6: Plot of factor heritabilities
# subplot(f1_row,f1_col,6)
# plot(G_h2s,'Color','r')
# hold on
# plot(E_h2s,'Color','b')
# hold off
# ylim([0 1])
# 
# 
# if sp_num>10,
# 
# %Figure of trace plots of the largest elements of each column of the factor loading matrix
# fig2 = figure(2);
# set(fig2,'Position',[0 800,600,800]);
# clf;
# f2_row=8;
# f2_col=2;
# for k=0:(min(2*f2_row,size(Posterior.Lambda,1)/p)-1)
# subplot(f2_row,f2_col,k+1)
# [~,o] = sort(-abs(mean(Posterior.Lambda((1:p)+k*p,max(1,sp_num-100):sp_num),2)));
# plot(Posterior.Lambda(o(1:5)+k*p,1:sp_num)')
#      end
#      
#      %Figure of trace plots and histograms of the factor heritablities
#      fig4 = figure(4);
#      set(fig4,'Position',[0 800,600,800]);
#      clf
#      f4_row=8;
#      f4_col=4;
#      for k=1:min(2*f4_row,size(Posterior.G_h2,1))
#      h2s = Posterior.G_h2(k,1:sp_num);
#      if sum(h2s(~isnan(h2s)))==0
#      continue
#      end
#      subplot(f4_row,f4_col,2*k-1)
#      plot(h2s)
#      %         ylim([0 max(sds(~isnan(sds)))])
#      ylim([0 1])
#      subplot(f4_row,f4_col,2*k)
#      hist(h2s,100)
#      xlim([0 1])
#      end               
#      
#      figure(fig1);  
#      
#      k = size(Posterior.Lambda,1)/p;
#      h2s = Posterior.G_h2(:,1:sp_num);     
#      G_Lambdas = zeros(size(Posterior.Lambda));
#      G_est = zeros(p,p);
#      E_est = zeros(p,p);
#      for j=1:sp_num
#      Lj = bsxfun(@times,reshape(Posterior.Lambda(:,j),p,k),1./sqrt(VY'));
# h2j = Posterior.G_h2(:,j);
# G_Lj = Lj * diag(sqrt(h2j));
# G_Lambdas(:,j) = G_Lj(:);
# Gj = G_Lj * G_Lj' + diag(1./(VY'.*Posterior.ps(:,j)));%
# G_est = G_est + Gj./sp_num;
# 
# E_Lj = Lj * diag(1-h2j)* Lj' + diag(1./(VY'.*Posterior.resid_ps(:,j)));
# E_est = E_est + E_Lj./sp_num;
# end
# G_Lambda = reshape(mean(G_Lambdas,2),p,k);
# 
# %plot 2: visualize posterior mean genetic correlation matrix
# subplot(f1_row,f1_col,2)
# clims=[-1 1];
# imagesc(CovToCor(G_est),clims)
# 
# %plot 4: visualize posterior mean genetic factor loading matrix
# subplot(f1_row,f1_col,4)
# clim=min(.75,max(max(abs(G_Lambda))));
# clims=[-clim clim];
# imagesc(G_Lambda(p_view,:),clims)
# lx = xlim;
# colorbar;  
# 
# %plot 6: posterior mean factor heritabilities
# subplot(f1_row,f1_col,6)
# plot(mean(h2s,2),'Color','r');
# ylim([0 1])
# hold on
# plot(1-mean(h2s,2),'Color','b');
# hold off
# colorbar
# xlim(lx)
# 
# drawnow()
# saveas(fig2,'Diagnostics_autoCorr_Lambda.jpg','jpg')
# saveas(fig4,'Diagnostics_factor_selection.jpg','jpg')
# end
# drawnow()
# saveas(fig1,'Diagnostics_Lambda.jpg','jpg')
# 
# end
# 
