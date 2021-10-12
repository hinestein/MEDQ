#Function that calculates the pth multivariate empirical dynamic quantile
#Inputs are X.list, p and method
#X.list is a list of nxT matrices for n time-series at T time intervals where each list entry is represents another dimension
#p is a point or a vector of quantiles to calculate, default is p = 0.5
#Method is determines which method of depth used. Can be Mahalanobis, Tukey, Liu or Oja.
#For dimensions greater than 2, Mahalanobis is by far the fastest
#Output is a point or vector row numbers for the pth quantile
#Can work for one dimensional time-series data
#Data set for each variable can be scaled in the function

MEDQ = function(X.list,p = 0.5, method = "Mahalanobis", scale = FALSE, weight = FALSE){
  d = length(X.list)
  if(!inherits(X.list,"list")){
    stop("X.list must be in the form of a list")
  }
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
  if(scale){
    for(i in 1:d){
      X.list[[i]] = scale(X.list[[i]])
    }
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
    for(j in 1:nrow(depth3)){
      z1 = x1[zt,d]
      if(d > 1){
        for(k in 1:(d-1)){
          z1 = z1 - (vec[k,1]/vec[d,1]) * (x1[j,k] - x1[zt,k])
        }
      }else{
        z1 = z1
      }
      if(x1[j,d] < z1){
        depth.dash[j,i] = depth3[j,i]
      }else{
        depth.dash[j,i] = 2 * depth3[zt,i] - depth3[j,i]
      }
    }
  }
  dist.1 = list()
  for(i in 1:nc){
    x1 <- NULL
    for(j in 1:d){
      x1 = cbind(x1, X.list[[j]][,i])
    }
    dist.1[[i]] = as.matrix(dist(x1, upper = TRUE, diag = TRUE))
  }
  n.3 = NULL
  pb <- txtProgressBar(min = 1, max = nr * nc, style = 3)
  t1 = 1
  n.2 = matrix(0, nrow = nr, ncol = length(p))
  for(i in 1:nr){
    n.1 = 0
    for(j in 1:nc){
      setTxtProgressBar(pb, t1)
      t1 = t1 + 1
      n.1 = n.1 + p * (dist.1[[j]][,i] %*% (depth.dash[,j] > depth.dash[i,j]))[1,1] + (1 - p) * (dist.1[[j]][,i] %*% (1 - (depth.dash[,j] > depth.dash[i,j])))[1,1]
    }
    n.2[i,] = n.1
  }
  n.3 = apply(n.2, 2, order)[1,]
  n.3
}

#Example for 2 dimensional positively correlated data
A.r1 = NULL
A.r2 = NULL
for(i in 1:500){
  rn = rnorm(1)
  A.r1 = rbind(A.r1, arima.sim(n = 200, list(ar = c(0.8897)),
                               sd = sqrt(sqrt(0.1796))) + rn)
  A.r2 = rbind(A.r2, arima.sim(n = 200, list(ar = c(0.8897)),
                               sd = sqrt(sqrt(0.1796))) + 2 * rn + 10)
}

X.list = list()
X.list[[1]] = A.r1
X.list[[2]] = A.r2

p = c(0.05,0.5,0.95)
M = MEDQ(X.list = X.list, p = p, method = "Mahalanobis")

par(mfrow = c(1,2))
ts.plot(t(X.list[[1]]))
lines(X.list[[1]][M[1],], col = "green")
lines(X.list[[1]][M[2],], col = "red")
lines(X.list[[1]][M[3],], col = "blue")

ts.plot(t(X.list[[2]]))
lines(X.list[[2]][M[1],], col = "green")
lines(X.list[[2]][M[2],], col = "red")
lines(X.list[[2]][M[3],], col = "blue")

