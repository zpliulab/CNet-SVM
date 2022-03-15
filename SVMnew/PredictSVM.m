function [y0, pred_class, AUC, RefereceResult] = PredictSVM(train_data, train_y, Est_theta) 
%%
% Est_theta = theta;
% train_data = X_test;
% train_y = y_test;
[n1,p1] = size(train_data); 
X = [ones(n1,1), train_data]; 
y0 = train_y; 
theta = Est_theta;
%% 2 different data possible: for GENES and for CLASSES (here only for classes)
% separating line: x*w+b=f(x)
pp = X*theta;
% AUC = plotroc(y0,pp);    % 对于SVM分类器，MATLAB有自己的自带方法plotroc方法，但是对于随机森林得到的分类模型和预测不适用
% AUC = AUC(y0,pp);
% pred_class = sign(pp);
pred_class = sign(pp);
figure;
AUC = plotROC(y0,pred_class);
%% 3. sensitivity and specificity  for CLASSES
tabClass = confusionmat(y0,pred_class);
figure;
tabPlot = printConMat(tabClass);
% stats = statsOfMeasure(tabClass);   % 有错误
[Result,RefereceResult] = confusion.getMatrix(y0,pred_class);
% disp(Result)
% disp(RefereceResult)
return