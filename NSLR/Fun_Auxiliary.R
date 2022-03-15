# 图正则项 \beta' M \beta 中的矩阵 M -------------------------------------------------------------------------
get.penalityMatrix=function(adj,X1, y1){
  feature.num = dim(adj)[2]       # M 为方阵，行数=列数
  M.c = diag(0,feature.num)       # M = 0
  M.lasso = diag(0,feature.num)   # M = 0
  M.elastic = diag(1,feature.num) # M = I 单位阵
  M.network = Non.NormalizedLaplacianMatrix(adj)      # M = Laplaceian 矩阵  
  M.AdaptNet = AdaptNet.Non.NormalizedLap(adj,X1, y1) # M = 自适应 Laplaceian 矩阵  
  return(list(M.c=M.c,M.lasso=M.lasso,M.elastic=M.elastic,M.network=M.network,M.AdaptNet=M.AdaptNet))
}
# 构造 network 矩阵-------------------------------------------------------------------------
# Normalized Laplacian Matrix from adjacency matrix
laplacianMatrix = function(adj){
  diag(adj) <- 0                   # 邻接矩阵对角元0
  # 度矩阵元素（对角）--邻接矩阵每行元素的绝对值之和 
  deg <- apply(abs(adj),1,sum)     # abs(adj)-矩阵各元素去绝对值、1-表示按行计算，2表示按列、sum-自定义的调用函数
  p <- ncol(adj)
  L <- matrix(0,p,p)               # p*p 的0元素的 Laplaceian 矩阵
  nonzero <- which(deg!=0)         # 哪些行 元素绝对值之和不为0
  for (i in nonzero){
    for (j in nonzero){
      L[i,j] <- -adj[i,j]/sqrt(deg[i]*deg[j])  # i j 不等时（L 为对称阵）
    }
  }
  diag(L) <- 1                                 # 对角线为1
  return(L)
}


# # 测试 ----------------------------------------------------------------------
# 
# X <- matrix(-1:-4,ncol=2) 
# X
# abs(X)
# adj <- X
# deg <- apply(abs(adj),1,sum)     
# nonzero <- which(deg!=0) 


# Non-Normalized Laplacian Matrix from adjacency matrix
Non.NormalizedLaplacianMatrix = function(adj){
  diag(adj) <- 0
  deg <- apply(adj,1,sum)
  D = diag(deg)
  L = D - adj             # 最普通的 L 矩阵 
  return(L)
}
# 求 AdaptNet 矩阵中的  正则的 M = L_star ------------------------------------------------------
AdaptNet.penality.matrix = function(adj, X, y){
  
  # ridge 的系数
  library(glmnet)
  glmnet.fit = glmnet(X, y, lambda=0, family='binomial')
  Beta = coef(glmnet.fit)
  
  p <- ncol(adj)
  coeff.sign = sign(Beta)[2:(p+1)]    # beta 的第一个元素为 截距，不是系数
  
  diag(adj) <- 0
  deg <- apply(abs(adj),1,sum)
  L <- matrix(0,p,p)
  nonzero <- which(deg!=0)
  for (i in nonzero){
    for (j in nonzero){
      temp.sign = coeff.sign[i]*coeff.sign[j]
      L[i,j] <- -temp.sign*adj[i,j]/sqrt(deg[i]*deg[j])
    }
  }
  diag(L) <- 1
  return(L_star=L)
}

# 求 AdaptNet 矩阵中的  非正则的 M = L_star
AdaptNet.Non.NormalizedLap = function(adj,X, y){
  library(glmnet)
  glmnet.fit = glmnet(X, y, lambda=0, family='binomial')
  Beta = coef(glmnet.fit)
  
  p <- ncol(adj)
  coeff.sign = sign(Beta)[2:(p+1)]
  diag(adj) <- 0
  deg <- apply(abs(adj),1,sum)
  L <- matrix(0,p,p)
  nonzero <- which(deg!=0)
  for (i in nonzero){
    for (j in nonzero){
      temp.sign = coeff.sign[i]*coeff.sign[j]
      L[i,j] <- -temp.sign*adj[i,j]
    }
  }
  diag(L) <- deg
  return(L_star=L)
}


# # 测试 ----------------------------------------------------------------------
# 
# Beta <- c(0, -1, 2)  # 0 为截距
# p <- 2
# coeff.sign = sign(Beta)[2:(p+1)]

# 计算 Sn、Se-------------------------------------------------------------------------
sen.spe = function(pred, truth){
  #----------------------------------------
  # 1: True sample
  # 0: False sample
  # truth = abs(sign(sim.data$w))
  # pred = abs(sign(out1$w[-length(tmp)]))
  #----------------------------------------
  sen = length(which(pred[which(truth==1)]==1))/length(which(truth==1))
  spe = length(which(pred[which(truth==0)]==0))/length(which(truth==0))
  result = c(sen,spe)
  names(result) = c("sensitivity","specificity")
  return(result)
}


# Prior network regularization
Prior.network = function(adj){
  # adj is a p x p matrix
  # L   is a (p+1) x (p+1) matrix
  L = laplacianMatrix(adj)
  L = rbind(L,rep(0,p))
  L = cbind(L,rep(0,p+1))
  return(L)
}

# -------------------------------------------------------------------------
get_sim_prior_Net = function(n,t,p11,p12){
  A = matrix(0,nrow=n,ncol=n) 
  for(i in 1:n){
    for(j in 1:n){
      if(i>j){
        set.seed(10*i+8*j)
        if(i<t&j<t){
          if(runif(1) < p11) A[i,j] = 1}    # runif(1)函数生成1个服从正态分布随机数
        else{
          if(runif(1) < p12) A[i,j] = 1}  
      }
    }
  }
  A = A + t(A)
  diag(A)=0
  return(A)
}

# -------------------------------------------------------------------------
GET.SIM.DATA = function(smaple.num, feature.num, random.seed, snrlam=0.05){
  # smaple.num = 700; feature.num = 100; random.seed = 10; snrlam=0.05
  ii = random.seed
  set.seed(30)
  w  <-c(rnorm(40),rep(0,(feature.num-40)))  # rnorm(40)产生40个服从正态分布的随机数
  b = 0
  mu <- rep(0,40)
  Sigma <- matrix(.6, nrow=40, ncol=40) + diag(40)*.4
  
  set.seed(ii*2)
  X1 <- mvrnorm(n=smaple.num, mu=mu, Sigma=Sigma)  # mvrnorm(n, mean, sigma)产生多元正态分布的随机数
  
  set.seed(ii*3)
  X2 <- matrix(rnorm(smaple.num*(feature.num-40), mean = 0, sd = 1), nrow = smaple.num, ncol = feature.num-40)
  X = cbind(X1,X2)

  Xw <- -X%*%w 
  pp <- 1/(1+exp(Xw))
  y  <- rep(1,smaple.num)
  y[pp<0.5] <- 0

  p = dim(X)[1]
  q = dim(X)[2]
  
  set.seed(ii*3)
  XX = X + snrlam*matrix(rnorm(p*q),ncol=q)
  
  return(list(X=XX,y=y,w=w))
}

# GET.SIM.DATA2 和 GET.SIM.DATA 没有区别-------------------------------------------------------------------------
# GET.SIM.DATA2 是2021.10.2，修改的，feature = 1000,增加非0变量feature.ture --------------------

GET.SIM.DATA2 = function(smaple.num, feature.num, feature.ture, random.seed, snrlam=0.05){
  # smaple.num = 700; feature.num = 100; random.seed = 10; snrlam=0.05
  ii = 10
  set.seed(30)
  w  <- c(rnorm(feature.ture),rep(0,(feature.num-feature.ture)));
  b = 0
  mu <- rep(0,feature.ture)
  Sigma <- matrix(.6, nrow = feature.ture, ncol = feature.ture) + diag(feature.ture)*.4
  
  set.seed(ii*2)
  X1 <- mvrnorm(n = smaple.num, mu = mu, Sigma = Sigma)
  
  set.seed(ii*3)
  X2 <- matrix(rnorm(smaple.num*(feature.num-feature.ture), mean = 0, sd = 1), nrow = smaple.num, ncol = feature.num-feature.ture)
  X = cbind(X1,X2)

  Xw <- -X%*%w 
  pp <- 1/(1+exp(Xw))
  y  <- rep(1,smaple.num)
  y[pp<0.5] <- 0
  
  p = dim(X)[1]
  q = dim(X)[2]
  
  set.seed(random.seed)
  XX = X + snrlam*matrix(rnorm(p*q),ncol=q)
  
  return(list(X=XX,y=y,w=w))
}