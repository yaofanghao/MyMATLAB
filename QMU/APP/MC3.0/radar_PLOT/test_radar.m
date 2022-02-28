% 绘制雷达图

% data = [5 8 1 4];
% lim = [0,10;0,10;0,10;0,10];
% labels = {'1','2','3','4'};
% Draw_radar(data,lim,labels);
% hold on;
% data = [7,5,4,8];
% lim = [0,10;0,10;0,10;0,10];
% labels = {'1','2','3','4'};
% Draw_radar(data,lim,labels);

data = [1.2 1.4 0 0 0];
lim = [0,3;0,3;0,3;0,3;0,3];
labels = {'1','2','3','4','5'};
legendlables = {'test','threshold'};
Draw_radar(data,lim,labels);
% hold on;
% data = [1,1,0,0,0];
% lim = [0,3;0,3;0,3;0,3;0,3];
% labels = {'1','2','3','4','5'};
% Draw_radar(data,lim,labels);