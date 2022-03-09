%--------------------------------------------------------------------------
%  examp13.4-1  读取examp13_4_1.xls中数据，进行因子分析
%--------------------------------------------------------------------------
% CopyRight：xiezhh

%*********************************读取数据*********************************
[X,textdata] = xlsread('examp13_4_1.xls');
X = X(:,3:end);
varname = textdata(4,3:end);
obsname = textdata(5:end,2);


%******************调用factoran函数根据原始观测数据作因子分析*****************
% 从原始数据（实质还是相关系数矩阵）出发，进行因子分析，公共因子数为4
% 进行因子旋转(最大方差旋转法)
[lambda,psi,T,stats] = factoran(X,4)
Contribut = 100*sum(lambda.^2)/8
CumCont = cumsum(Contribut)



% 从原始数据（实质还是相关系数矩阵）出发，进行因子分析，公共因子数为2
% 进行因子旋转(最大方差旋转法)
[lambda,psi,T,stats,F] = factoran(X, 2)
Contribut = 100*sum(lambda.^2)/8
CumCont = cumsum(Contribut)

[varname' num2cell(lambda)]

%**************将因子得分F分别按耐力因子得分和速度因子得分进行排序*************
obsF = [obsname, num2cell(F)];
F1 = sortrows(obsF, 2);    % 按耐力因子得分排序
F2 = sortrows(obsF, 3);    % 按速度因子得分排序
head = {'国家/地区','耐力因子','速度因子'};
result1 = [head; F1];
result2 = [head; F2];

%*************************绘制因子得分负值的散点图***************************
plot(-F(:,1),-F(:,2),'k.');
xlabel('耐力因子得分（负值）');
ylabel('速度因子得分（负值）');
box off;
gname(obsname);
