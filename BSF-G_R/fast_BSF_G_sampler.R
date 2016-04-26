require(matlab)
require(PEIP)
require(Matrix)
require(pracma)

fast_BSF_G_sampler = function(burn, sp, thin, b0, b1, h2_divisions, epsilon, priors, draw_iter){
  #% -- Daniel Runcie -- %

  #% Gibbs sampler for genetic covariance estimation based on mixed effects
  #% model, with missing data
  
  #% Based on:
   
  #% Runcie and Mukherjee (2013) Dissecting high-dimensional traits
  #% with Bayesian sparse factor analysis of genetic covariance matrices.
  #% GENETICS.

  #% (c) April 22, 2013

  #% code modified from original provided by Anirban Bhattacharya
  
  #% This function implements the BSF-G partially collapsed Gibbs sampler. It loads all input data
  #% and matrices from setup.mat in the current directory. Priors and control
  #% parameters are passed to the function.

  #% setup.mat is a struct with at least:
  #%     Y: data matrix
  #%     X: fixed effect design matrix
  #%     Z_1: random effect design matrix for factor model
  #%     Z_2: additional random effect design matrix
  #%     A: Additive genetic relationship matrix
   
  #%    For analysis of Ayroles et al 2009 data, can include:
  #%     Ayroles_results: struct holding gene names and correlations estimated in that paper
     
  #%    For analysis of simulations:
  #%     U_act: r x p matrix of true breeding values
  #%     E_act: n x p matrix of true model residuals
  #%     gen_factor_Lambda: p x k_G matrix of true genetic factor loadings
  #%     error_factor_Lambda: p x k matrix of true residual factor loadings
  #%     h2: p x 1 vector of true heritabilities
  #%     factor_h2s: k x 1 vector of true latent factor heritabilities
  #%     G, R: p x p matrix of true genetic and residual covariances
 
  #% The function takes the following inputs:
  #%     burn: number of burnin samples
  #%     sp: total number of samples to collect
  #%     thin: thinning rate of chain
  #%     b0,b1: parameters controlling rate of adaptation of factor model size
  #%     h2_divisions: number of discrete steps for each factor heritability parameter
  #%     epsilon: truncation point for factor loadings during adaptation
  #%     draw_iter: frequency of updating diagnostic plots
  #%     priors: struct holding various prior hyperparameters: 
  #%         k_init: initial number of factors to initialize
  #%         as, bs: inverse gamma hyperparameters for model residuals, as well as non-factor random effects
  #%         df: degrees of freedom for t-distribution of ARD prior on factor loadings
  #%         ad1, bd1: inverse gamma hyperparamters for first factor shrinkage multiplier (/delta_1)
  #%         ad2, bd2: inverse gamma hyperparamters for remaining factor shrinkage multiplier (/delta_i, i \in 2...k)
 
  #% The function output is the struct, Posterior, with the following fields:
  #  %     Lambda, U, no_f, ps, resid_ps, delta, G_h2: matrices with each column a 
  #%         (vectorized if necessary) posterior sample of the:
  #  %             Lambda: factor loading matrix
  #%             U: genetic random effect matrix
  #%             no_f: number of significant factors
  #%             ps: genetic residual precision after accounting for latent factors
  #%             resid_ps: phenotypic residual precision
  #%             delta: column shrinkage parameters for Lambda
  #%             G_h2: factor heritabilties
  #%     B, d, W: matrices of posterior means of fixed effect (B), 
  #%         residual genetic (d), or 2nd random effect (W) coefficients. Some may be empty
  #% 
  #% Several diagnostic plots are produced during the run. 
  #%     Their interpretation is described within the source codes:
  #%         draw_simulation_diagnostics.m: For simulated data with known true values
  #%         draw_results_diagnostics.m: Otherwise

  clear Y
  clear Z_1
  clear Z_2
  clear X
  
  global Y   #%n x p matrix of phenotypic data
  global Z_1   #%n x r incidence matrix for additive genetic effects
  global Z_2   #%n x r2 incidence matrix for another set of random effects
  global X   #%n x b design matrix of fixed effects
  
  
  nrun = burn + sp * thin #% number of posterior samples
  k_min = 1e-1 #% minimum factor loading size to report in running status
  prop = 1.00  #% proportion of redundant elements within columns necessary to drop column

  #% ------read data--------%
  load('../setup.mat')

  #%Determine if 'setup.mat' contains output of a simulation, based on if
  #%known factor loadings are included. Affects plotting functions
  simulation = true
  if(isempty(objects('gen_factor_Lambda'))){
    simulation = false
  }
  simulation

  #%normalize Y to have zero mean and unit variances among observed values,
  #%allowing for NaNs.
  c(n,p) = dim(Y)
  Y_full = Y
  Mean_Y = zeros(1, ncol(Y))
  VY = zeros(1, ncol(Y))
  
  for(j in 1:p){
    Mean_Y[j] = mean(Y[~is.nan(Y[,j]),j])
    VY[j] = var(Y[~is.nan(Y[,j]),j])
    if(is.nan(VY[j])){
      VY[j] = 1
    }
  }
  
  Y = Y - Mean_y                     
  Y = Y * (1/sqrt(VY))      

  #%determine if a design matrix (X) exists (is loaded from setup.mat). If
  #%not, make a dummy X-matrix with no columns.
  if(~exists('X','var')){
    X = zeros(0,n)
  }
  if(ncol(X) ~= n){
    X = zeros(0,n)
  }

  #%Determine if a second random effects design matrix exists. If not, make a
  #%dummy matrix
  if(~exists('Z_2','var')){ #exists(Z_2) & exists('var')
    Z_2 = zeros(0,n)
  }
  if(ncol(Z_2) ~= n){
    Z_2 = zeros(0,n)
  }

  #%calculate model dimensions
  r = nrow(Z_1)
  r2 = nrow(Z_2)
  b = nrow(X)


  #% --- Initialize variables --- %
  #  %residual parameters. This structure holds the priors hyperparamters for
  #%the gamma prior on the model residual variances. It also holds the current
  #%estimate of the residual precision
  resid$as = priors$as
  resid$bs = priors$bs
  resid$Y = Y
  resid$p = p
  resid$ps = rgamma(p, rate = resid$as, scale = 1/resid$bs)           #%residual precision

  #%Factors. This struct holds all information about the latent factors,
  #%including the current number of factors, priors hyperparameters for the
  #%factor Loadings, as well as current values of the factor loadings, their
  #%precision, the factor scores, and the genetic heritability of each factor
  k = priors$k_init      #% initial number of factors
  df = priors$df           #% prior degrees of freedom for hierarchical t-prior on loadings 
  ad1 = priors$ad1         #% priors on delta
  bd1 = priors$bd1
  ad2 = priors$ad2
  bd2 = priors$bd2
  Factors$r = n            
  Factors$n = n
  Factors$p = p
  Factors$k = k
  Factors$df = df
  Factors$ad1 = ad1
  Factors$bd1 = bd1
  Factors$ad2 = ad2
  Factors$bd2 = bd2
  Factors$psijh = matrix(rgamma(p*k, rate = df/2, scale = 2/df), nrow = p, ncol = k)   #%individual loadings precisions
  Factors$delta = [gamrnd(ad1+10,1/bd1);gamrnd(ad2,1/bd2,[k-1,1])];  #% components of tauh !!!!!!need size
  Factors$tauh = cumprod(Factors$delta)                              #%extra shrinkage of each loading column
  Factors$Plam = Factors$psijh * t(Factors$tauh)          #%total precision of each loading
  Factors$Lambda = zeros(p,k) + matrix(runif(m*n), nrow = m, ncol = n) * reshape(sqrt(1/Factors$Plam),p,k)   #%factor loadings
  Factors$h2 = runif(k)                                             #%factor heritability
  Factors$h2_divisions = h2_divisions                                #%discretizations of heritability
                        
  Factors$num = 0
  Factors$no_f = zeros(sp,1)
  Factors$nofout = k * ones(nrun,1)
                        
  #%genetic_effects. This structure holds information about latent genetic
  #%effects. U is latent genetic effects on factor traits. d is genetic
  #%effects on residuals of the factor traits. Plus prior hyperparameters for
  #%genetic effect precisions
                        genetic_effects$n = nrow(Z_1)
                        as = priors$as
                        bs = priors$bs
                        genetic_effects$as = as
                        genetic_effects$bs = bs
                        genetic_effects$ps = rgamma(p, rate = as, scale = 1/bs)
                        genetic_effects$U = matrix(rnorm(k*r),k,r) * sqrt(Factors$h2)
                        genetic_effects$d = matrix(rnorm(p*r),p,r) * 1/sqrt(genetic_effects$ps)
                        
                        #%interaction_effects. Similar to genetic_effects structure except for
                        #%additional random effects that do not contribute to variation in the
                        #%factor traits
                        as = priors$as
                        bs = priors$bs
                        interaction_effects$as = as
                        interaction_effects$bs = bs
                        interaction_effects$ps = rgamma(p, rate = as, scale = 1/bs)
                        interaction_effects$mean = zeros(p,r2)
                        interaction_effects$n = r2
                        interaction_effects$W = matrix(rnorm(p*r2),p,r2) * 1/sqrt(interaction_effects$ps)
                        interaction_effects$W_out = zeros(p,r2)
                        
                        #%fixed_effects hold B
                        fixed_effects$b = b
                        fixed_effects$cov = zeros(b,b) #%inverse covariance of fixed effects
                        fixed_effects$mean = zeros(p,b) #%mean of fixed effects
                        fixed_effects$B = matrix(rnorm(p*b),p,b) #%current estimate of fixed effects
                        
                        Factors$scores = genetic_effects$U * Z_1 + matrix(rnorm(k*n),k,n) * sqrt(1-Factors$h2) #%initialize factor scores
                        
                        #%Posterior holds Posterior samples and Posterior means
                        Posterior$Lambda = zeros(0,sp)
                        Posterior$no_f = zeros(sp,1)
                        Posterior$ps = zeros(p,sp)
                        Posterior$resid_ps = zeros(nrow(resid$ps),sp)
                        Posterior$B = zeros(p,nrow(X))
                        Posterior$U = zeros(0,sp)
                        Posterior$d = zeros(p,nrow(Z_1))
                        Posterior$W = zeros(p,nrow(Z_2))
                        Posterior$delta = zeros(0,sp)
                        Posterior$G_h2 = zeros(0,sp)
                        
                        #%save run parameters and hyperparameters
                        
                        params$p = p
                        params$n = n
                        params$r = r
                        params$Mean_Y = Mean_Y
                        params$VY = VY
                        params$b0 = b0
                        params$b1 = b1
                        params$epsilon = epsilon
                        params$prop = prop
                        params$as = priors$as
                        params$bs = priors$bs
                        params$df = priors$df
                        params$ad1 = priors$ad1
                        params$bd1 = priors$bd1
                        params$ad2 = priors$ad2
                        params$bd2 = priors$bd2
                        params$burn = burn
                        params$thin = thin
                        params$sp = sp
                        params$nrun = nrun
                        params$h2_divisions = h2_divisions
                        
                        if(simulation){
                          params$U_act = U_act
                          params$Lambda = error_factor_Lambda
                          params$h2 = h2
                          params$G = G
                          params$R = R
                          params$B = B
                          params$factor_h2s = factor_h2s
                          params$name = name
                        }
                        
                        #%precalculate some matrices
                        #%invert the random effect covariance matrices
                        Ainv = solve(A)
                        A_2_inv = diag(1, nrow(Z_2)) #%Z_2 random effects are assumed to have covariance proportional to the identity. Can be modified.
                        
                        #%pre-calculate transformation parameters to diagonalize aI + bZAZ for fast
                        #%inversion: inv(aI + bZAZ) = 1/b*u*diag(1./(s+a/b))*u'
                        #%uses singular value decomposition of ZAZ for stability when ZAZ is low
                        #%rank
                        ZAZ = t(Z_1) %*% A %*% Z_1
                        [u,s,~] = svd(ZAZ)$d
                        eig_ZAZ$vectors = u
                        eig_ZAZ$values = diag(s)
    
                        #%fixed effects + random effects 1
                        #%diagonalize mixed model equations for fast inversion: 
                        #%inv(a*blkdiag(fixed_effects.cov,Ainv) + b*[X; Z_1][X; Z_1]') = Q*diag(1./(a.*s1+b.*s2))*Q'
  Design = rbind(X, Z_1)
  Design2 = Design %*% t(Design)
  [~,~,q,S1,S2] = GSVD(Cholesky(blkdiag(fixed_effects$cov, Ainv)), Cholesky(Design2))
  svd_Design_Ainv$Q = t(solve(q))
  svd_Design_Ainv$s1 = diag(t(S1) %*% S1)
  svd_Design_Ainv$s2 = diag(t(S2) %*% S2)
  Qt_Design = t(svd_Design_Ainv$Q) %*% Design      
                          
                          #%random effects 2
                          #%as above, but for random effects 2. Here, fixed effects will be conditioned on, not sampled simultaneously. Otherwise identical.
                          Design = Z_2
                          Design2 = Design %*% t(Design)
                          [~,~,q,S1,S2] = GSVD(Cholesky(A_2_inv), Cholesky(Design2))
                          svd_Z2_2_A2inv$Q = t(solve(q))
                          svd_Z2_2_A2inv$s1 = diag(t(S1) %*% S1)
                          svd_Z2_2_A2inv$s2 = diag(t(S2) %*%S2)
                          Qt_Z2 = t(svd_Z2_2_A2inv$Q) %*% Design
                          
                          #%genetic effect variances of factor traits
                          #% diagonalizing a*Z_1*Z_1' + b*Ainv for fast inversion
                          #%diagonalize mixed model equations for fast inversion: 
                          #% inv(a*Z_1*Z_1' + b*Ainv) = Q*diag(1./(a.*s1+b.*s2))*Q'
                          #%similar to fixed effects + random effects 1 above, but no fixed effects.
                                  ZZt = Z_1 %*% t(Z_1)
                                  [~,~,q,S1,S2] = GSVD(Cholesky(ZZt), Cholesky(Ainv))
                                  svd_ZZ_Ainv$Q = t(solve(q))
                                  svd_ZZ_Ainv$s1 = diag(t(S1) %*% S1)
                                  svd_ZZ_Ainv$s2 = diag(t(S2) %*% S2)
                                  
                                  #%------start gibbs sampling-----%
                                  sp_num = 0
                                  tic
                                  for(i in 1:nrun){
                                    #%fill in missing phenotypes
                                    #%conditioning on everything else
                                    phenMissing = is.nan(Y_full) # which indices
                                    if(sum(sum(phenMissing)) > 0){
                                    meanTraits = fixed_effects$B %*% X +  genetic_effects$d %*% Z_1 
                                    + interaction_effects$W %*% Z_2 + Factors$Lambda %*% Factors$scores
                                    meanTraits = t(meanTraits)        
                                    resids = matrix(randn(dim(Y_full)), nrow(Y_full), ncol(Y_full)) * 1/sqrt(t(resid$ps))
                                    Y[phenMissing] = meanTraits[phenMissing] + resids[phenMissing]
                                    }

  #%sample Lambda
  #%conditioning on W, X, F, marginalizing over D
  Ytil = t(Y) - fixed_effects$B %*% X - interaction_effects$W %*% Z_2
  [Factors] = sample_lambda(Ytil, Factors, resid, genetic_effects, eig_ZAZ)
  
  #%sample fixed effects + random effects 1 ([B;D])
  #%conditioning on W, F, L
  Ytil = t(Y) - interaction_effects$W %*% Z_2 - Factors$Lambda %*% Factors$scores
  N = genetic_effects$n + fixed_effects$b
  [location_sample] = sample_means(Ytil, Qt_Design, N, resid, genetic_effects$ps, svd_Design_Ainv) #CHECK
  fixed_effects$B = location_sample[,1:fixed_effects$b]
  genetic_effects$d = location_sample[,fixed_effects$b+1:fixed_effects$b+genetic_effects$n]
  
  #%sample random effects 2
  #%conditioning on B, D, F, L
  Ytil = t(Y) - fixed_effects$B %*% X - genetic_effects$d %*% Z_1 - Factors$Lambda %*% Factors$scores
  N = interaction_effects$n
  if(N > 0){
    location_sample = sample_means(Ytil, Qt_Z2, N, resid, interaction_effects$ps, svd_Z2_2_A2inv)
    interaction_effects$W = location_sample
  }

  #%sample factor h2
  #%conditioning on F, marginalizing over U
  Factors = sample_h2s_discrete(Factors, eig_ZAZ)
  
  #%sample genetic effects (U)
  #%conditioning on F, Factor h2
  genetic_effects = sample_Us(Factors, genetic_effects, svd_ZZ_Ainv)
  
  #%sample F
  #%conditioning on U, Lambda, B, D, W, factor h2s
  Ytil = t(Y) - fixed_effects$B %*% X - genetic_effects$d %*% Z_1 - interaction_effects$W %*% Z_2
  Factors = sample_factors_scores(Ytil, Factors, resid, genetic_effects)
  
  #% -- Update ps -- %
  Lambda2 = (Factors$Lambda)^2    
  Factors$psijh = matrix(rgamma(nrow(Factors$df)*ncol(Factors$df), shape = Factors$df/2 + 0.5, scale = 2/(Factors$df + Lambda2 * t(Factors$tauh))), nrow(Factors$df), ncol(Factors$df))
                                                                  
  #%continue from previous Y residual above
  Ytil = Ytil - Factors$Lambda %*% Factors$scores
  n = nrow(Y)
  resid$ps = matrix(rgamma(nrow(resid$as)*ncol(resid$as), shape = resid$as + 0.5*n, scale = 1/(resid$bs+0.5 %*% t(rowSums(Ytil^2)))), nrow(resid$as), ncol(resid$as))  #%model residual precision
  n = genetic_effects$n
  genetic_effects$ps =  matrix(rgamma(nrow(genetic_effects$as)*ncol(genetic_effects$as), shape = genetic_effects$as + 0.5*n, scale = 1/(genetic_effects$bs+ 0.5 %*% t(rowSums(genetic_effects$d^2)))),
                               nrow(genetic_effects$as), ncol(genetic_effects$as)) #%random effect 1 (D) residual precision
  n = interaction_effects$n
  interaction_effects$ps = matrix(rgamma(nrow(interaction_effects$as)*ncol(interaction_effects$as), shape = interaction_effects$as + 0.5*n, scale = 1/(interaction_effects$bs+0.5 %*% t(rowSums(interaction_effects$W^2)))),
                                  nrow(interaction_effects$as), ncol(interaction_effects$as)) #%random effect 2 (W) residual precision
                                                                    
  #%------Update delta & tauh------%
  c(delta, tauh) = sample_delta(Factors, Lambda2)
  Factors$delta = delta
  Factors$tauh = tauh
                                                                  
  #%---update precision parameters----%
  Factors$Plam = Factors$psijh * t(Factors$tauh)
                                                
  #% ----- adapt number of factors to samples ----%
  c(Factors, genetic_effects) = update_k(Factors, genetic_effects, b0, b1, i, epsilon, prop)
                                              
  #% -- save sampled values (after thinning) -- %
  if((i %% thin) == 0 && i > burn){
                                                
    sp_num = (i-burn)/thin
                                                  
    Posterior =  save_posterior_samples(sp_num, params, Posterior, 
                                        resid,fixed_effects, genetic_effects, Factors, 
                                        interaction_effects)
                                                
    if(mod(sp_num, 100) == 0){
      save('Posterior', 'Posterior', 'params')
    }
  }
                                              
  #% -- provide run diagnostics and plots -- %
  if((i %% draw_iter) == 0){   
    directory = strread(pwd, '%s', 'delimiter', '/')
    disp(directory(end))
    disp(i)
    Factors$nofout[i] - Factors$num
    elapsed = toc
  
    #%output some running statistics on the current factors and their
    #%genetic variances
    c(Factors$delta, t(c(1:Factors$k)), Factors$h2, t(sum(t(Factors$scores)^2))/(nrow(Y)-1), t(sum(t(genetic_effects$U)^2))/(nrow(Z_1)-1))
    disp(strcat('Time remaining:', num2str((nrun-i) * (elapsed/i) * 1/60)))
  
    #%make some plots of some running statistics
    if(simulation){
      draw_simulation_diagnostics(i,sp_num,params,Factors,genetic_effects,resid,Posterior,gen_factor_Lambda,error_factor_Lambda,G,R,h2)
    }
    else{
      draw_results_diagnostics(i, sp_num, params, Factors, Posterior)
    }
  }
}
  toc
  save('Posterior','Posterior','params')
  
  return(Posterior, params)
}