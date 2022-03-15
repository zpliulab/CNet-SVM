%% 
clc; clear; close all; format long;  
%% 
global L1  X_train  y_train  lambda_opt alpha_opt m q
%% 
adj = importdata('adj_example.csv');
A_adj = sparse(adj.data);
%% Set folds in cv
nfold = 5;
%% 
vector_hat = importdata('vector_hat_example.txt'); 
delta_hat = vector_hat.data;
%% Set options
options = optimoptions('fmincon','Algorithm','interior-point');
options.MaxFunEvals = 1e7;
options.MaxIterations = 1e4;
%% 
txt = importdata(['train10.txt']);  
train_data = txt.data;
[m,q] = size(train_data);
X_train = train_data(1:m,2:q);    
% y_train = train_data(1:m,1);    % 0- -1,  1 - 1;
y_train = 2*(train_data(1:m,1)+1)-3;    % 0- -1,  1 - 1;
p = q-1;
%% laplace 
[L, L1] = Laplacian_Matrix(p,A_adj);
%% candidate lambdas in simulations
% e=(log(100)-log(.0001))/99;
% lambda=exp(log(.0001):e:log(100)); 
%% set candidate values of alpha 
% alpha = 0.1:.2:.9;  
%% provide the value of alpha
alpha = 0.5;   
%% get the least upper bound, lambda_max
lammax = getLambMaxSVM(X_train, y_train, alpha); 
e = (log(lammax)-log(1))/19;
lambda = exp(log(1):e:log(lammax)); 
%% get the optimal lambda and alpha through cross validation 
[lambda_opt, alpha_opt, r] = cvSVM( X_train, y_train, L1, alpha, lambda, nfold );
%% optimal parameters
alpha_opt 
lambda_opt  
 
%% Restrictions
theta_0 = zeros(q,1);  
u_0 = (1e-5)*ones(q,1);  

a11 = eye(q,q);
a22 = -1*eye(q,q);
A1 = [a11;a22];
A2 = [a22;a22];
A = [A1,A2];
b = zeros(2*q,1);

Aeq = []; 
beq = [];
  
vlb1 = zeros(q,1);        
vub1 = zeros(q,1);  
delta1 = [0; abs(delta_hat)];
for j = 1:q
    if delta1(j) == 1
        vlb1(j) = 1e-3;
        vub1(j) = 1;
    else 
        vlb1(j) = -inf;
        vub1(j) =  inf;
    end
end
vlb = [vlb1; vlb1];        
vub = [vub1; vub1];  
%% interior point method
[X_sol, cost, exitflag, output, mu, grad, hessian] = fmincon(@(x)(costFunctionSVM(x(1:q),x(q+1:2*q))), [theta_0;u_0], A, b, Aeq, beq, vlb, vub, [], options);
theta = X_sol(1:q);
u = X_sol(q+1: 2*q);
theta_hat = X_sol(2:q)  
theta_hat0 = X_sol(1);   
u_hat0 = X_sol(q+1);
%% print
fprintf('Cost at theta found by fmincon: %e', cost);
fprintf('\n')
%% predict
format long;
txt1 = importdata(['test500.txt']);  
test_data = txt1.data;
[k,l] = size(test_data);
X_test = test_data(1:k,2:l);   
y_test = 2*(test_data(1:k,1)+1)-3; 
[y_true, y_pre, AUC, RefereceResult] = PredictSVM(X_test, y_test, theta); 
disp(RefereceResult)
P_test = [y_pre, y_true];
para = [lambda_opt, AUC, cost]
%% save
csvwrite(['theta_hat_500.csv'],theta_hat);
csvwrite(['theta_500.csv'],theta);
csvwrite(['P_test_500.csv'],P_test);
csvwrite(['para_500.csv'],para);

