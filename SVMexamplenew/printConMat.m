% https://download.csdn.net/download/miao_5/10201008?spm=1001.2101.3001.5697
function plot = printConMat(conMat)
%load result.mat conMat;

[m,n] = size(conMat);
x = 1:m;
y = 1:n;

imagesc(conMat);        %# Create a colored plot of the matrix values 

%% 我自己加的 2021.8.28  https://zhuanlan.zhihu.com/p/365519245
% label = {'-1','1'};
xlabel('Predict class','fontsize',12);
ylabel('Actual class','fontsize',12);
%%
colormap(flipud(gray)); %# Change the colormap to gray (so higher values are  
                        %#   black and lower values are white)  
  
textStrings = num2str(conMat(:),'%0.02f');      %# Create strings from the matrix values  
textStrings = strtrim(cellstr(textStrings));    %# Remove any space padding  

[sx,sy] = meshgrid(x,y);                      %# Create x and y coordinates for the strings  

hStrings = text(sx(:),sy(:),textStrings(),...     %# Plot the strings  
                'HorizontalAlignment','center');  
midValue = mean(get(gca,'CLim'));               %# Get the middle value of the color range  
textColors = repmat(conMat(:) > midValue,1,3);  %# Choose white or black for the  
                                                %#   text color of the strings so  
                                                %#   they can be easily seen over  
                                                %#   the background color  
set(hStrings,{'Color'},num2cell(textColors,2)); %# Change the text colors  
%% 自己修改的横纵坐标的标签
% plot = set(gca,...
%         'XTick',1:m,...                         %# Change the axes tick marks  
%         'XTickLabel',label,...  %#   and tick labels  
%         'YTick',1:n,...  
%         'YTickLabel',label,...  
%         'TickLength',[0 0]);  

%%
xStrings = strcat('c',num2str(x(:),'%-2d'));
yStrings = strcat('c',num2str(y(:),'%-2d'));

plot = set(gca,...
        'XTick',1:m,...                         %# Change the axes tick marks  
        'XTickLabel',xStrings,...  %#   and tick labels  
        'YTick',1:n,...  
        'YTickLabel',yStrings,...  
        'TickLength',[0 0]);  
%%
set(gca,'FontSize',12,'LineWidth',1.5);
set(get(gca,'YLabel'),'FontSize',12);
set(get(gca,'XLabel'),'FontSize',12);

end