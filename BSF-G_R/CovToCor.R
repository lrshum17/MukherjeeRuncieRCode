CovToCor = function(X){
  corX = diag(1/sqrt(diag(X))) %*% X %*% diag(1/sqrt(diag(X)))
  return(corX)
}