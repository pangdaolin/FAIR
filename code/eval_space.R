#Vector and trace correlation coefficients (Ye and Weiss, JASA, 2003)
eval.space <- function(A, B, orthnm = TRUE) 
{
  if(!is.matrix(A)) A <- as.matrix(A)
  if(!is.matrix(B)) B <- as.matrix(B)
  
  if(orthnm)
  { 
    A <- qr.Q(qr(A))
    B <- qr.Q(qr(B)) 
  }
  
  mat <- t(B) %*% A %*% t(A) %*% B
  d <- eigen(mat)$values
  d <- (d + abs(d))/2
  q <- sqrt(prod(d))
  r <- sqrt(mean(d))
  
  c(q, r, acos(q))
}

