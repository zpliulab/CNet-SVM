% predict-分类器对测试集的分类结果
% ground. _truth -测试集的正确标签,这里只考虑二分类,即-1和1
% auc-返回ROC曲线的曲线 下的面积
function auc = plotROC(ground_truth, predict)
%初始点为( 1.0, 1.0 )
%计算出ground_ truth中正样本的数目pos_ num和负样本的数目neg. num
pos_num = sum(ground_truth==1);
neg_num = sum(ground_truth==-1);
m=size(ground_truth,1);
[pre,Index]= sort(predict);
ground_truth=ground_truth(Index);
x=zeros(m+1,1);
y=zeros(m+1,1);
auc=0;
x(1)=1;
y(1)=1; 
for i=2:m
    TP=sum(ground_truth(i:m) == 1);
    FP= sum(ground_truth(i:m) == -1);
    x(i)=FP/neg_num;
    y(i)=TP/pos_num;
    auc=auc + (y(i)+y(i-1))*(x(i-1)-x(i))/2;
end
x(m+1) = 0;
y(m+1) = 0;
auc = auc+y(m)*x(m)/2;
% plot(x,y);

h1 = plot(x,y,'-b','LineWidth',3,'MarkerSize',3);
%% 自己加的
hold on;
color4 = [107,105,102]./255;
x_dig=0:0.01:1;
y_dig=x_dig;
plot(x_dig,y_dig,'--','Color',color4,'LineWidth',1.5); %画中间的虚线

xlabel('False Positive Ratio (1-Specificity)','fontsize',12,'FontWeight','bold');
ylabel('True Positive Ratio (Sensitivity)','fontsize',12,'FontWeight','bold');
% title('ROC曲线图');

set(gca,'FontSize',12,'LineWidth',1.5);
set(get(gca,'YLabel'),'FontSize',12);
set(get(gca,'XLabel'),'FontSize',12);
set(gca,'YTick',[0:0.2:1]);

grid on
ROCtitle_1=['AUC = ',num2str(roundn(auc,-3))];
hh = legend([h1], ROCtitle_1, 'Location','southeast')%,'Location','southeast');
set(hh,'edgecolor','white');
end

%%
% function  auc = plot_roc( predict, ground_truth )
% % INPUTS
% %  predict       - 分类器对测试集的分类结果
% %  ground_truth - 测试集的正确标签,这里只考虑二分类，即0和1
% % OUTPUTS
% %  auc            - 返回ROC曲线的曲线下的面积
% 
% %初始点为（1.0, 1.0）
% x = 1.0;
% y = 1.0;
% %计算出ground_truth中正样本的数目pos_num和负样本的数目neg_num
% pos_num = sum(ground_truth==1);
% neg_num = sum(ground_truth==0);
% %根据该数目可以计算出沿x轴或者y轴的步长
% x_step = 1.0/neg_num;
% y_step = 1.0/pos_num;
% %首先对predict中的分类器输出值按照从小到大排列
% [predict,index] = sort(predict);
% ground_truth = ground_truth(index);
% %对predict中的每个样本分别判断他们是FP或者是TP
% %遍历ground_truth的元素，
% %若ground_truth[i]=1,则TP减少了1，往y轴方向下降y_step
% %若ground_truth[i]=0,则FP减少了1，往x轴方向下降x_step
% for i=1:length(ground_truth)
%     if ground_truth(i) == 1
%         y = y - y_step;
%     else
%         x = x - x_step;
%     end
%     X(i)=x;
%     Y(i)=y;
% end
% %画出图像     
% plot(X,Y,'-ro','LineWidth',2,'MarkerSize',3);
% xlabel('虚报概率');
% ylabel('击中概率');
% title('ROC曲线图');
% %计算小矩形的面积,返回auc
% auc = -trapz(X,Y);  
% 
% end

% https://blog.csdn.net/xmu_jupiter/article/details/21885299