clc; clear; close all; format long;
%% 初始化
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

%% 图例
% set(h00,'handlevisibility','off') % 只要想要的图例
% set(h000,'handlevisibility','off') % 只要想要的图例
legend('Lasso-SVM','ENet-SVM (\alpha=0.5)','SCAD-SVM (a=3.7)', 'L2SCAD-SVM (a=3.7, \alpha=0.5)',...
        'location', 'north'); % eastoutside,在图外
xlabel('\theta');
ylabel('P(\theta) ');
% title('L_{q} penalties');

%%    0 -尽量不与数据冲突，自动放置在最佳位置 ; 1 -放置在图形的右上角 ; 2 -放置在图形的左上角 ;  
%%    3 -放置在图形的左下角 ;4 -放置在图形的右下角 ; -1 -放置在图形视窗的外右边 ;
%%    lgd=legend('第一个图例','第二个图例','orientation','horizontal'); %默认为vertical，改为horizontal
%% 图中图
% magnify


    
