MEDQ = function(X.list,p = 0.5, method = "Mahalanobis"){
  d = length(X.list)
  nr = NULL
  for(i in 1:d){
    nr = c(nr, nrow(X.list[[i]]))
  }
  if(length(unique(nr)) > 1){
    stop("Dimensions are not equal")
  }
  nr = nr[1]
  nc = NULL
  for(i in 1:d){
    nc = c(nc, ncol(X.list[[i]]))
  }
  nc = nc[1]
  if(length(unique(nc)) > 1){
    stop("Dimension are not equal")
  }
  if(max(p) > 1 | min(p) < 0){
    stop("p must be between 0 and 1")
  }
  if(method == "Mahalanobis"){
    depth3 = NULL
    for(i in 1:nc){
      my_sample <- NULL
      for(j in 1:d){
        my_sample = cbind(my_sample, X.list[[j]][,i])
      }
      mu <- apply(my_sample,2,mean)
      sigma <- cov(my_sample)
      m_dist <- mahalanobis(my_sample, mu, sigma)
      m_depth <- 1/(1 + m_dist)
      depth3 = cbind(depth3, m_depth)
    }
  }else{
    depth3 = matrix(0, nrow = nr, ncol = nc)
    for(i in 1:nc){
      my_sample <- NULL
      for(j in 1:d){
        my_sample = cbind(my_sample, X.list[[j]][,i])
      }
      for(j in 1:nr){
        depth3[j,i] = depth(my_sample[j,], my_sample, method = method)
      }
    }
  }
  depth.dash = matrix(0, nrow = nr, ncol = nc)
  for(i in 1:nc){
    x1 <- NULL
    for(j in 1:d){
      x1 = cbind(x1, X.list[[j]][,i])
    }
    zt = which(depth3[,i] == max(depth3[,i]))[1]
    sigma <- cov(x1)
    EVV=eigen(sigma)
    vec=EVV$vectors
    lambda=EVV$values
    for(j in 1:nrow(depth3)){
      z1 = x1[zt,d]
      for(k in 1:(d-1)){
        z1 = z1 - (vec[k,1]/vec[d,1]) * (x1[j,k] - x1[zt,k])
      }
      if(x1[j,d] < z1){
        depth.dash[j,i] = depth3[j,i]
      }else{
        depth.dash[j,i] = 2 * depth3[zt,i] - depth3[j,i]
      }
    }
  }
  dist.1 = list()
  for(i in 1:ncol(A.r1)){
    x1 <- NULL
    for(j in 1:d){
      x1 = cbind(x1, X.list[[j]][,i])
    }
    dist.1[[i]] = as.matrix(dist(x1, upper = TRUE, diag = TRUE))
  }
  n.3 = NULL
  for(k in p){
    n.2 = NULL
    for(i in 1:nrow(A.r1)){
      n.1 = 0
      for(j in 1:ncol(A.r1)){
        n.1 = n.1 + k * dist.1[[j]][,i] %*% (depth.dash[,j] > depth.dash[i,j]) + (1 - k) * dist.1[[j]][,i] %*% (1 - (depth.dash[,j] > depth.dash[i,j]))
      }
      n.2 = c(n.2, n.1)
    }
    n.3 = c(n.3, which(n.2 == min(n.2)))
  }
  n.3
}