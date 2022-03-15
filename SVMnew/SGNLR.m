% SGNLR = function(X, y, M, lambda0, alpha, niter=100, gate = 0){
%   # X = X1  
%   # y = y1  
%   # 
%   # M = PM$M.network
%   # # M = PM$M.lasso
%   # lambda0 = 0.2 
%   # alpha  = 0.3
%   # 
%   # # M = PM$M.c 
%   # # lambda0 = 0  
%   # # alpha = 0  
%   # niter=20  
%   # gate = 0
%   
%   # -----------------
%   # X   is a n x (p+1) matrix
%   # y   is a n x 1 vector
%   # M   is a (p+1) x (p+1) matrix
%   # -----------------
%   if(!is.matrix(y)){y = matrix(y,ncol = 1)}
%   if(!is.matrix(X)){X = matrix(y,ncol = ncol(X))}
%   
%   n = dim(X)[1]
%   p = dim(X)[2]
%   
%   # Adjusting penalty parameters
%   lambda = n*alpha*lambda0
%   eta = n*(1-alpha)*lambda0
%   
%   X = cbind(X,rep(1,n))
%   # X   is a n x (p+1) matrix
%   
%   if(dim(M)[1]!= p+1){
%     M = cbind(M,rep(0,p)); 
%     M = rbind(M,rep(0,p+1))
%     } 
%   # M   is a (p+1) x (p+1) matrix
%   # -----------------
%   # Initialization
%   w = matrix(0, nrow = p+1)
%   err=0.0001
%   Z = NULL
%   
%   obj=log.likelihhod = NULL
%   for(i in 1:niter){
%     # i = 1
%     w0 = w
%     D = update.diagMatrix(X, w, gate)
%     z = update.zVector(X, y, w, D)
%     XDX = t(X)%*%D%*%X + eta*M
%     t   = t(X)%*%D%*%z
%     for(j in 1:p){
%       # j = 1
%       w[j] = 0
%       w[j] = (Soft.thresholding(t[j]-t(w)%*%XDX[j,], lambda))/XDX[j,j]
%       
%       if (w[j] != 0){Z = c(Z,j)}
%     }
%     
%     w[p+1] = 0
%     w[p+1] = (t[p+1]-t(w)%*%XDX[p+1,])/XDX[p+1,p+1]
%     
%     # Calculate J (for testing convergence)
%     log.likelihhod = c(log.likelihhod,objective(X,y,w))
%     obj = c(obj,logreg.obj(X,y,w,M,lambda0, alpha))
%     
%     if(i>3&&abs(obj[i]-obj[i-1])/abs(obj[i-1])<0.00001){break}
%   }
%   names(w) = c(paste("w", 1:(length(w)-1), sep = ""), "Intercept")
%   return(list(w=w, log.likelihhod = log.likelihhod, obj=obj, Z=Z))
% }


% function [theta, Z] = SGNLR(X, y, M, lambda0, alpha, niter, gate) % [theta, log_likelihhod, obj, Z] 
% 
X = train_data; 
y = train_y;  
M = M_network;
lambda0 = 0.2;
alpha  = 0.3;
niter=20;
gate = 0;

[n,p] = size(X);

% Adjusting penalty parameters
lambda = n*alpha*lambda0;
eta = n*(1-alpha)*lambda0;

X = [X,ones(n,1)];
[n1, p1] = size(M);
if n1~= p+1
    M = [M,zeros(p,1)];   % 左右合并
    M = [M;zeros(p+1,1)']; % 上下合并
end
M_hat = M;    % 原M为sparse，不可行
for t = 1:p
    M_hat(t,t) = 0;  
end
M_hat = abs(M_hat);

% h1=view(biograph(M_hat));   %画图
% Z = graphtraverse(M_hat,1); %使用深度优先算法从第1个节点开始遍历
% order2=graphtraverse(M_hat,4,'Method','BFS'); %使用广度优先遍历


% Initialization
theta = zeros(p+1,1);
err=0.0001;

Z = zeros(p, 1);
obj = zeros(niter, 1);;
log_likelihhod = zeros(niter, 1);;


 for i = 1:niter
% i = 1;
    theta0 = theta;
    D = Update_Diag_Matrix(X, theta, gate);
    z = Update_z_Vector(X, y, theta, D);
    XDX = X'*D*X + eta*M;
    t   = X'*D*z;
    for j = 1:p
%       j = 1;
      theta(j) = 0;
      D_hat = XDX;
      theta(j) = (Soft_Thresholding(t(j)-theta'*D_hat(j,:)', lambda))/XDX(j,j);
%       [ci sizes] = components(M); 
     
% j = 1;      % 初始时，从第一行开始
% path = [j];
%     for k=1:p   
%         %寻找当前行k中的1所在的位置     
%         for l=1:p       
%             if(M_hat(j,l) == 1)            
%                 j = l;            
%                 path = [path, j];            
%                 break;        
%             end
%         end
%     end
    
      if M(i,i) ~=0
          for k = 2:p
              if M(i,k) ~= 0;
                 Z(k) = k;
              else
                  break % 必须加 break
              end
          end
      end
%       if theta(j) ~= 0
%           Z(j,2) = j;
%       end
%       if M(i,i) ~=0 && theta(j) ~= 0
%           Z(j,3) = j;
%       end      
    end
    theta(p+1) = 0;
    theta(p+1) = (t(p+1)-theta'*D_hat(p+1,:)')/D_hat(p+1,p+1);
    
    % Calculate J (for testing convergence)
    log_likelihhod(i) = Objective(X,y,theta);
    obj(i) = Logreg_Objective(X,y,theta,M,lambda0,alpha);
    if i>3 && abs(obj(i)-obj(i-1))/abs(obj(i-1)) < 0.00001
        break
    end
 end
       
% path
% return

