function [theta_hat, theta_0] = LogitisLapSVM (X,y,L,lambda1,lambda2)
% global p 
%% To get the theta which minimize the cost function -sum_i{1-y_i(theta'*x_i)}
%% +lambda1 |theta|_1+lambda2 theta' A theta
%% The input y is a n by 1 matrix
%% Input X is G X n gene expression matrix
%% L is graphical Laplacian matrix, lambda1 and lambda2 are parameters in the cost function that are known
%% The output theta is the optimal solution of the covex optimization problem
%% The dimension of theta: G by 1
% X = X_train;
% y = y_train;

[dimn,~]=size(y); % get the number of n and set it to dimn
[dimG,m]=size(X); % get the number of G and set it to dimG

% size(ones(dimn,1))
% size(X_train)
% aa = [ones(dimn,1), X];


%%  To use cvx to solve the convex problem
cvx_begin
    variables theta(dimG+1);
%   minimize (  sum((1-y)'.*(theta'*X))+sum( log_sum_exp( [zeros(1,m); -theta'*X]) ) +lambda1*norm(theta,1)+lambda2*theta'*L*theta  );
%   minimize (  sum((1-y)'.*(theta'*X))+sum( log_sum_exp( -theta'*X ) ) +lambda1*norm(theta,1)+lambda2*theta'*L*theta  );
%     sum( max( zeros(dimn,1), ones(dimn,1)-y.*( [ones(dimn,1), X] * theta) ) ) + lambda1*norm(theta(2:dimG+1),1)+lambda2*theta(2:dimG+1)'*L*theta(2:dimG+1);
%     sum( max( zeros(dimn,1), ones(dimn,1)-y.*( [ones(dimn,1), X] * theta) ) );
%     minimize( sum( max( zeros(dimn,1), ones(dimn,1)-y.*( [ones(dimn,1), X] * theta) ) ) );
%% 测试 penalty 函数，可行
%     minimize( lambda1*norm(theta(2:dimG+1),1)+lambda2*theta(2:dimG+1)'*L*theta(2:dimG+1) );
%% 测试 max 函数，可行
    m = dimn;
    y_train = y;
    X_trainhat = [ones(m,1), X'];   % 这是因为cvSVM line 68 是转置，所以这里X'
    z = X_trainhat * theta;
%     minimize( sum( max( zeros(m,1), ones(m,1)-y_train.*z ) ) );   % 这是把costFunction 13 的目标函数拿过来了
%% 在加上惩罚项
    minimize( sum( max( zeros(m,1), ones(m,1)-y_train.*z ) ) + lambda1*norm(theta(2:dimG+1),1)+lambda2*theta(2:dimG+1)'*L*theta(2:dimG+1) );
%% 先固定，lambda/alpha, 运行interior 算法，得到theta，进行测试，为了解决 “horzcat 串联的矩阵的维度不一致”。
% X = X_train;
% y = y_train;
% size(y);
% size( [ones(dimn,1), X] * theta);
% aaa = [ones(dimn,1), X] * theta;
% bbb = y.*( [ones(dimn,1), X] * theta);
% ccc = ones(dimn,1)-y.*( [ones(dimn,1), X] * theta);
% ddd = max(zeros(dimn,1), ccc);
% sum(ddd);
% size(y.*( [ones(dimn,1), X] * theta));
% sum( max( zeros(dimn,1), ones(dimn,1)-y.*( [ones(dimn,1), X] * theta) ) )
%%

cvx_end
theta_hat = theta(2:dimG+1);
theta_0 = theta(1);

return

    
    
