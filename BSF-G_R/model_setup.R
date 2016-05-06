## prepare output directory
rm(list = ls())
a = as.double(format(Sys.time(), "%OS1")) #reset random seeds
s = 1e6*(a %% 1)

set.seed(s)
rm(Posterior) #clear previous runs

rep = 1 #can be set by shell at runtime for repeated runs.
folder = paste0('rep', rep)
dir.create(folder)
setwd(folder)

#if needed, add path to model scripts and functions
#addpath(genpath('../../scripts'))




## Initialize run parameters and prior hyper-parameters
# run parameters 

draw_iter=100
burn=1000
sp=1000
thin=10
b0=1
b1=0.0005
epsilon=1e-2
h2_divisions = 100

# prior hyperparamters
k_init = 20 # initial number of factors to initialize
as=2 # inverse gamma hyperparameters for model residuals, as well as non-factor random effects
bs=1/10 # inverse gamma hyperparameters for model residuals, as well as non-factor random effects
df = 3 # degrees of freedom for t-distribution of ARD prior on factor loadings
ad1 = 2.1 # inverse gamma hyperparamters for first factor shrinkage multiplier (/delta_1)
bd1 = 1/20 # inverse gamma hyperparamters for first factor shrinkage multiplier (/delta_1)
ad2 = 3 # inverse gamma hyperparamters for remaining factor shrinkage multiplier (/delta_i, i \in 2...k)
bd2 = 1 # inverse gamma hyperparamters for remaining factor shrinkage multiplier (/delta_i, i \in 2...k)

save(as,bs, k_init,df,ad1,bd1,ad2,bd2,file='priors.Rdata')
priors = load('priors.Rdata')




## Run Gibbs sampler
source("fast_BSF_G_sampler.R")
Posterior = fast_BSF_G_sampler(burn,sp,thin,b0,b1,h2_divisions,epsilon,priors,draw_iter)


## Calculate posterior means of G, R, P, Lambda, factor_h2s, h2s
p = nrow(Posterior$ps)
sp_num = ncol(Posterior$Lambda)
VY = matrix(1,1,p)
k = nrow(Posterior$Lambda)/p
h2s = Posterior$G_h2[,1:sp_num]  
G_Lambdas = matrix(0, nrow(Posterior$Lambda), ncol(Posterior$Lambda))

Lambda_est = matrix(0,p,k)
G_est = matrix(0,p,p)
P_est = matrix(0,p,p)
factor_h2s_est = matrix(0,k,1)

for (j in 1:sp_num){
  Lj = matrix(Posterior$Lambda[,j],p,k) * (1/sqrt(t(VY)))
  h2j = Posterior$G_h2[,j]
  
  Pj = Lj %*% t(Lj) + diag(1/(t(VY) * Posterior$ps[,j])) + diag(1/(t(VY) * Posterior$resid_ps[,j]))
  Gj = Lj %*% diag(h2j) %*% t(Lj)  + diag(1/(t(VY) * Posterior$ps[,j]))
  
  Lambda_est = Lambda_est + Lj/sp_num
  P_est = P_est  + Pj/sp_num
  G_est = G_est + Gj/sp_num
  
  factor_h2s_est = factor_h2s_est + h2j/sp_num
}

source("cor.R")
cors = abs(cor(t(Lambda_est),t(params$Lambda)))

k = ncol(params.Lambda)
o = matrix(0,k,1)
for (j in 1:k){
  o[j] = which.max(cors[,j])
}
# o=order(max_cors)
Lambda_est_ordered = Lambda_est[,o]
factor_h2s_est = factor_h2s_est[o]

factor_angles = matrix(0,k,1)
for (j in 1:k){
  a = Lambda_est_ordered[,j]
  b = params.Lambda[,j]
  costheta = sum(a*b)/(sqrt(sum(a^2))*sqrt(sum(b^2)))
  theta = acos(costheta)
  factor_angles[j] = pi/2-abs(pi/2-theta)
}

# Figure of trace plots and histograms of the factor heritablities
f4_row=8;
f4_col=4;
par(mfrow=c(f4_row,f4_col))
for (ji in 1:min(4*8/2,length(o))){
  j = o(ji)
  h2s = Posterior$G_h2[j,]
  if (sum(h2s[!isnan(h2s)])==0){
    next
  }
  plot(h2s, ty="l", ylim = c(0,1))
  hist(h2s, 100, xlim= c(-0.1,1))
}

posterior_mean$G = G_est
posterior_mean$P = P_est
posterior_mean$Lambda = Lambda_est
posterior_mean$factor_order = o
posterior_mean$factor_h2s_est = factor_h2s_est
posterior_mean$factor_angles = factor_angles

save('posterior_mean',file='Posterior_mean.Rdata')

## Produce plot of results
# figure(5)
p = nrow(Posterior$ps); 
sp_num = ncol(Posterior$Lambda)
VY = matrix(1,1,p)
k = nrow(Posterior$Lambda)/p
h2s = Posterior$G_h2[,1:sp_num]    
G_Lambdas = matrix(0, nrow(Posterior$Lambda), ncol(Posterior$Lambda))
G_est = matrix(0,p,p)
R_est = matrix(0,p,p)
for (j in 1:sp_num){
  Lj = matrix(Posterior$Lambda[,j],p,k) * (1/sqrt(t(VY)))
  h2j = Posterior$G_h2[,j]
  G_Lj = Lj %*% diag(sqrt(h2j))
  G_Lambdas[,j] = matrix(G_Lj, nrow=nrow(G_Lj)*ncol(G_Lj) ncol=1)
  Gj = G_Lj %*% t(G_Lj) + diag(1/(t(VY) * Posterior$ps[,j]))
  G_est = G_est + Gj/sp_num
  
  E_Lj = Lj %*% diag(1-h2j) %*% t(Lj) + diag(1/(t(VY) * Posterior$resid_ps[,j]))
  R_est = R_est + E_Lj/sp_num
}

G_Lambda = matrix(apply(G_Lambdas,1,mean),p,k)

f1_row = 1
f1_col=2
par(mfrow=c(f1_row,f1_col))
clims=c(-1, 1)
image(CovToCor(G_est),zlim=clims)

clim=min(.75,max(abs(G_Lambda)))
clims=c(-clim, clim)
image(G_Lambda,zlim=clims)
# Check lx = xlim;
# Check colorbar; 

# export_fig('transparent','append','nocrop','/Users/der7/Documents/Statistics/Sparse_factor_G_matrix/Genetics_paper/BGSFM/Instructions/Sim.pdf')

setwd('..')