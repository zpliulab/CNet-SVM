clc; clear; close all; format long;
%% ��ʼ��
N = 10;
l = 4;

% N = 20;
% l = 1;
n = 1/N;
m = 2;
x = -l:n:l;

%% L1
y1 = abs(x);
h1 = plot(x, y1, '--xr', 'linewidth', 2.5);
hold on

%% L1 + L2
y12 = 0.5 .* x .^2 + 0.5 * abs(x);
h12 = plot(x, y12, '--xm','linewidth', 2.5);
hold on;

%% SCAD
lambda = 1;
a = 3.7;
yscad = lambda * abs(x) .* (abs(x) <= lambda) + ...
    + ( -(x.^2 - 2*a*lambda*abs(x) + lambda^2)/(2*(a-1)) ) .* (lambda < abs(x) & abs(x) <= a*lambda) + ...
    + (a+1)*lambda^2/2 .* (abs(x)> a*lambda);
h_scad = plot(x, yscad, 'g','linewidth', 2.5);
hold on

%% SCAD
lambda = 1;
a = 3.7;
yscad = 0.5 .* x .^2 + 0.5*(lambda * abs(x) .* (abs(x) <= lambda) + ...
    + ( -(x.^2 - 2*a*lambda*abs(x) + lambda^2)/(2*(a-1)) ) .* (lambda < abs(x) & abs(x) <= a*lambda) + ...
    + (a+1)*lambda^2/2 .* (abs(x)> a*lambda));
h_scad = plot(x, yscad, 'b','linewidth', 2.5);
hold off

%% ͼ��
% set(h00,'handlevisibility','off') % ֻҪ��Ҫ��ͼ��
% set(h000,'handlevisibility','off') % ֻҪ��Ҫ��ͼ��
legend('Lasso-SVM','ENet-SVM (\alpha=0.5)','SCAD-SVM (a=3.7)', 'L2SCAD-SVM (a=3.7, \alpha=0.5)',...
        'location', 'north'); % eastoutside,��ͼ��
xlabel('\theta');
ylabel('P(\theta) ');
% title('L_{q} penalties');

%%    0 -�����������ݳ�ͻ���Զ����������λ�� ; 1 -������ͼ�ε����Ͻ� ; 2 -������ͼ�ε����Ͻ� ;  
%%    3 -������ͼ�ε����½� ;4 -������ͼ�ε����½� ; -1 -������ͼ���Ӵ������ұ� ;
%%    lgd=legend('��һ��ͼ��','�ڶ���ͼ��','orientation','horizontal'); %Ĭ��Ϊvertical����Ϊhorizontal
%% ͼ��ͼ
% magnify


    
