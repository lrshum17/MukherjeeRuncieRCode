sample_delta = function(Factors, Lambda2_resid){
  #sample delta and tauh parameters that control the magnitudes of higher
  #index factor loadings.
  
  # Reading the input, its a structure element
  ad1 = Factors$ad1
  ad2 = Factors$ad2
  bd1 = Factors$bd1
  bd2 = Factors$bd2
  k = Factors$k
  delta = Factors$delta
  tauh = Factors$tauh
  psijh = Factors$psijh
  
  # Apply element-by-element binary operation to two arrays with singleton expansion enabled
  # In this case ints doing array multiplication of psijh and Lambda2_resid
  # mat = bsxfun(@times,psijh,Lambda2_resid);
  mat = psijh * Lambda2_resid
  
  #Finding the number of rows of the matrix mat
  n_genes = dim(mat)[1]
  
  # Doing some matrix addition and multiplication
  ad = ad1 + 0.5 * n_genes * k
  
  # delta(1) is the 1st element of the vector delta
  bd = bd1 + 0.5 * (1/delta[1]) * sum(tauh * sum(t(mat)))
  
  #  gamrnd is a random number generationg see website
  # http://www.mathworks.com/help/stats/gamrnd.html
  delta[1] = rgamma(1, shape = ad, scale = 1/bd)
  
  # Finds the relative product of the vector delta 
  # i.e if a vector A = [1 2 3 4 5] then
  # B = cumprod(A) = [1 1*2 1*2*3 1*2*3*4 1*2*3*4*5] = [1 2 6  24 120]
  # http://www.mathworks.com/help/matlab/ref/cumprod.html
  tauh = cumprod(delta)
  
  # For looping for values of h, 2 to k
  for(h in 2:k){
    #n_genes*(k-h+1) is selecting the h-k+1 index value form the vector n_genes
    ad = ad2 + 0.5 * n_genes[k-h+1] 
    
    # the function sum is a summation of the vector elements 
    # delta(h) is picking the h element from the vector delta
    # sum(tauh(h:end) is summing the values h to end index in tauh
    # the use of ' means its taking a transpose
    bd = bd2 + 0.5 * (1/delta[h]) %*% sum(tauh[h:length(tauh)] * sum(t(mat[,h:length(tauh)])))
    
    # Again generating a random number 
    delta[h] = rgamma(1, shape = ad, scale = 1/bd)
    
    tauh = cumprod(delta)
  }
  
  return(list(delta, tauh))
}