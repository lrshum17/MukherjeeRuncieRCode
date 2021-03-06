Cor = function(X, Y){
  #correlation matrix of correlations between rows of X and Y
  
  # Finding the number of rows of each matrix X and Y
  x = dim(X)[1]
  y = dim(Y)[1]
  
  # Concatenating matrices, stacking maxtirx X on top of matrix Y
  X = rbind(X, Y)
  
  # Finds the convultion/covariance of the tranpose of X
  # X' --> transpose of X
  covX = cov(t(X))
  
  # Initialzies a diagonal matrix with the elements of sqrt(1./diag(covX))
  # sqrt(1./diag(covX)) has to be a vector of elements, 
  # 1.diag(conX) is done element wise 
  # diag(conX) is extracting the diagonal elements of the matrix X
  d = diag(sqrt(1/diag(covX)))
  
  # Matrix multiplication
  r = d %*% covX %*% d
  
  # Extracts elements of the r matrix, to make a new matrix,
  # extracts the elements of rows 1 to x and in columns x+1 to x+y, 
  # recall that x and y the number of rows of the X and Y matrix
  r = r[1:x, (x+1):(x+y)]
  return(r)
}