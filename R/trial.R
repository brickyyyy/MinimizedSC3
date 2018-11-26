calculate_consensus<-fuction(matrix, k){
  #constructing a binary matrix for the cluster identities n
  n = ncol(matrix)
  c = nrow(matrix)
  b = matrix(0L,nrow = c,ncol = n)
  message("Calculating consensus matrix...")
  
  for (i in 1:n) {
    for (j in 1:c) {
      #populate every cell of the binary matrix
      rowValue=matrix[j,i]+k*(j-1)
      b[i,rowValue]=1
    }
  }
  
  inputMatrix=t(b)/n*c
  #add tolerance at convergence=1e-10.
  res=kmeans(x=inputMatrix, centers = k)
  return (res)
}