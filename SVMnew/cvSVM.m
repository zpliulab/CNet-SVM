function [lambda_opt, alpha_opt, r] = cvSVM(X,y,L,alpha,lambda,nfold)
%% The function for logistic regression network regularization applications  
%% The function is to use cross validation (CV) to get the optimal value of
%% parameter lambda and alpha through the training and testing procedures.
%% nfold is the number of folds in CV, which could be selected randomly
%% For Lasso case, alpha is 1, which is included in the code instruction.
% X = X_train;
% y = y_train;
% L = L1;

X = X';
[dimG,~]=size(X); % get the dimension G 
[dimn,~]=size(y); % get the dimension n
%% Get the number of candidate lambda's
[m,n]=size(lambda);
if m>n
    sizelam=m;
else
    sizelam=n;
end
%% Get the length of alpha. for lasso: alpha=1
[m,n]=size(alpha);
if m>n
    sizealpha=m;
else
    sizealpha=n;
end

%% To divide X and y into several folds respectively
 for i=1:(nfold-1)
      fnw=['x',int2str(i),'=X(:,(floor(dimn/nfold)*',int2str(i),'-floor(dimn/nfold)+1):(floor(dimn/nfold)*',int2str(i),'));'];
      eval(fnw);
      fnw=['y',int2str(i),'=y((floor(dimn/nfold)*',int2str(i),'-floor(dimn/nfold)+1):(floor(dimn/nfold)*',int2str(i),'),1);'];
      eval(fnw);
 end
 
 %% The last one is the columns that left 
 fnw=['x',int2str(nfold),'=X(:,(floor(dimn/nfold)*',int2str(nfold-1),'+1):dimn);'];
 eval(fnw);
 fnw=['y',int2str(nfold),'=y((floor(dimn/nfold)*',int2str(nfold-1),'+1):dimn);'];
 eval(fnw);
 
 %% To set the head and rear xnot first.
 xnot1=X(:,floor(dimn/nfold)+1:dimn);
 fnw=['xnot',int2str(nfold),'=[X(:,1:(floor(dimn/nfold)*',int2str(nfold-1),'))];'];
 eval(fnw);
 
 %% Set the rest xnot's
 for i=2:(nfold-1)
    fnw=['xnot',int2str(i),'=[X(:,1:(floor(dimn/nfold))*(',int2str(i),'-1))  X(:, (floor(dimn/nfold)*',int2str(i),'+1):dimn)];'];
    eval(fnw);
 end
 
  %% Set ynot in a similar way
 ynot1=y(floor(dimn/nfold)+1:dimn,1);
 fnw=['ynot',int2str(nfold),'=[y(1:(floor(dimn/nfold)*',int2str(nfold-1),'),1)];'];
 eval(fnw);
 for i=2:(nfold-1)
    fnw=['ynot',int2str(i),'=[y(1:(floor(dimn/nfold))*(',int2str(i),'-1)); y((floor(dimn/nfold)*i+1):dimn, 1)];'];
    eval(fnw);
 end
 
 %% The following is the cross validation to choose optimal alpha and lambda
for k=1:sizealpha
 alphacan = alpha(k);    
 for i=1:nfold
     for j=1:sizelam
         fnw=['[Theta,Theta0]=LogitisLapSVM (xnot',int2str(i),',ynot',int2str(i),', L,lambda(',int2str(j),')*alphacan,lambda(',int2str(j),')*(1-alphacan));'];
         eval(fnw);
         %% r is the matrix that contains the residual=|| yi-(1/(1+exp(-theta'*Xi))>0.5) ||_2
%       fnw=['r(',int2str(i),',',int2str(j),')=norm((transpose(y',int2str(i),')-(1./(1+exp(-transpose(Theta)*x',int2str(i),'))>.5)),2);'];
        fnw=['r(',int2str(i),',',int2str(j),')=norm((transpose(y',int2str(i),')-sign(transpose(Theta)*x',int2str(i),'+ Theta0)),2);'];
%        size(transpose(Theta))
%        size(x1)
%        size(y1)
        %%
        eval(fnw);
     end
 end
 
 %% to get the minimal residual and the corresponding index, which is the index for optimal lambda
 rmin(k)= max(mean(r,1)); % first set the maximal mean to rmax
 optindex=1;  % in case, the optimal lambda is the first one
 for j=1:sizelam
     if rmin>mean(r(:,j))
         rmin(k)=mean(r(:,j)); 
         optindex=j;  % if mean is less, take the index and reset the minimal residual
     end
 end
 %% Get the optimal lambda
 lambdaopt(k)=lambda(optindex);  % For every alpha, there will be an optimal lambda 
end
[err Ind]=min(rmin); 
alpha_opt=alpha(Ind);
lambda_opt=lambdaopt(Ind);

    