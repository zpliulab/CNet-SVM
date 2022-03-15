function lammax = getLambMaxSVM(X,y,alpha)
[n p] = size(X);
%% the first term in lambda_max computation
% a1 = 2*norm(X'*y ,inf);
a1 = norm(X'*y ,inf);
%% the second part of lambda_max calculation 
% Xi = zeros(1,p);
% for i = 1:n
%     Xi = Xi+X(i,:);
% end
% a2 = norm(Xi',inf);

%% to get lammax
% lammax = (a1+a2)/(2*alpha);
lammax = a1/alpha;
return
