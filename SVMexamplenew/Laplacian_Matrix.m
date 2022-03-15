% Non.NormalizedLaplacianMatrix = function(adj){
%   diag(adj) <- 0
%   deg <- apply(adj,1,sum)
%   D = diag(deg)
%   L = D - adj             # 最普通的 L 矩阵 
%   return(L)
% }

function [L, L1] = Laplacian_Matrix(p,adj) % Non Normalized

for i = 1:p
    adj(i,i) = 0;
end

deg = sum(adj,2);   % 矩阵各行元素求和
D = sparse(p,p);
for i = 1:p
    D(i,i) = deg(i);
%     D(i,i) = sqrt(deg(i));
end
L = D - adj;
L = sparse(L);

% D = diag(sum(adj,2)); 
D1 = sqrt(D); 
[dimd,~] = size(D1);
for i=1:dimd
    if ( D1(i,i)~= 0 )
        D1(i,i)=1/D1(i,i);
    end
end
L1 = D1*L*D1;
L1 = sparse(L1);
return
