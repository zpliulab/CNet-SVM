%% 格式化
clc; clear; close all; format long;  
%% 全局变量
global L1  X_train  y_train  lambda_opt alpha_opt m q
%% 邻接矩阵
adj = importdata('adj_example.csv');
A_adj = sparse(adj.data);
%% Set folds in cv
nfold = 5;
%% 割点
vector_hat = importdata('vector_hat_example.txt'); 
delta_hat = vector_hat.data;
%% Set options
options = optimoptions('fmincon','Algorithm','interior-point');
options.MaxFunEvals = 1e7;
options.MaxIterations = 1e5;
%% 循环 30 次
% for i = 1
% i = 12; 
%% 数据输入
txt = importdata(['train260.txt']);  
train_data = txt.data;
[m,q] = size(train_data);
X_train = train_data(1:m,2:q);     % 以索引的前1000个数据点作为测试样本Xtest
% y_train = train_data(1:m,1);     % 0- -1,  1 - 1;
y_train = 2*(train_data(1:m,1)+1)-3;    % 0- -1,  1 - 1;
p = q-1;
%% laplace 矩阵
[L, L1] = Laplacian_Matrix(p,A_adj);
% L = eye(p,p);
%% candidate lambdas in simulations
% e=(log(100)-log(.0001))/99;
% lambda=exp(log(.0001):e:log(100)); 
%% set candidate values of alpha 
% alpha = 0.2:.2:.6;  
%% provide the value of alpha
% alpha = 1;     % --- lasso
% alpha = 0.5;   % -- elastic net
%% get the least upper bound, lambda_max
% lammax = getLambMaxSVM(X_train, y_train, alpha); 
% e = (log(lammax)-log(1))/19;
% lambda = exp(log(1):e:log(lammax)); 
%% get the optimal lambda and alpha through cross validation 
% [lambda_opt, alpha_opt, r] = cvSVM( X_train, y_train, L1, alpha, lambda, nfold );
%% 最优参数
alpha_opt = 0.5
lambda_opt = 100
 
%% 约束条件
theta_0 = zeros(q,1); % 初值为0
u_0 = (1e-5)*ones(q,1);  % 初值为1

a11 = eye(q,q);
a22 = -1*eye(q,q);
A1 = [a11;a22];
A2 = [a22;a22];
A = [A1,A2];
b = zeros(2*q,1);

% A = [];
% b = [];
Aeq = []; 
beq = [];
% vlb = [];        
% vub = [];  

vlb1 = zeros(q,1);        
vub1 = zeros(q,1);  
delta1 = [0; abs(delta_hat)];
for j = 1:q
    if delta1(j) == 1
        vlb1(j) = 5e-5;
        vub1(j) = 1;
    else 
        vlb1(j) = -inf;
        vub1(j) =  inf;
    end
end
vlb = [vlb1; vlb1];        
vub = [vub1; vub1];  
%% 内点法求解
[X_sol, cost, exitflag, output, mu, grad, hessian] = fmincon(@(x)(costFunctionSVM(x(1:q),x(q+1:2*q))), [theta_0;u_0], A, b, Aeq, beq, vlb, vub, [], options);
theta = X_sol(1:q);
u = X_sol(q+1: 2*q);
theta_hat = X_sol(2:q)  % 系数
theta_hat0 = X_sol(1);  % 截距
u_hat0 = X_sol(q+1);
%% 输出 Cost 和 theta
fprintf('Cost at theta found by fmincon: %e', cost);
fprintf('\n')
%% 预测
format long;
txt1 = importdata(['test260.txt']);  
test_data = txt1.data;
[k,l] = size(test_data);
X_test = test_data(1:k,2:l);   
% y_test = test_data(1:k,1); 
y_test = 2*(test_data(1:k,1)+1)-3; 
[y_true, y_pre, AUC, RefereceResult] = PredictSVM(X_test, y_test, theta); 
disp(RefereceResult)
% [AUC_test, y_pre, y_true] = PredictSVM(X_train, y_train, theta);  
P_test = [y_pre, y_true];
para = [lambda_opt, AUC, cost]
%% 保存
csvwrite(['theta_hatopt55examplelambda10.csv'],theta_hat);
csvwrite(['thetaopt55examplelambda10.csv'],theta);
csvwrite(['P_testopt55examplelambda10.csv'],P_test);
csvwrite(['paraopt55examplelambda10.csv'],para);
% end
