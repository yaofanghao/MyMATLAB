clc
clear all
%%
%1.matlab数据统计分析常用函数
%（1）.计算样本均值
% mean(X)     % X为向量，则返回X的均值
% mean(A)     % A为矩阵，则返回每列的均值
% mean(A,2)   % A为矩阵，则返回每行的均
% 2）求样本方差
% var(X)       % X为向量，则返回X的方差
% var(A,[],1)  % A为矩阵，则返回每列的方差
% var(A,[],2)  % A为矩阵，则返回每行的方差
% 3）求标准差
% std(X)         % X为向量，则返回X的标准差
% std(A,flag,1)  % A为矩阵，则返回每列的标准差(flag取0或1，取0表示自由度=样本个数-1，取1时表示自由度=样本个数)
% std(A,flag,2)  % A为矩阵，则返回每行的标准差(flag取0或1，取0表示自由度=样本个数-1，取1时表示自由度=样本个数)
% %% （4）.求协方差
% cov(X)       % X为向量，则返回X的协方差
% var(A)       % A为矩阵，则返回该矩阵的协方差矩阵，其对角线元素为原矩阵A各列的方差
% cov(X,Y)     % 两等长列向量X、Y的协方差
%（5）.求相关系数
% corrcoef(X,Y)   % 等长列向量X、Y的相关系数
% corrcoef(A)     % 矩阵A的列向量的相关系数矩阵
%%
%2.不同概率模型抽样函数
% normrnd 可以生成一定均值和标准差的正态分布
% gamrnd 可以生成gamma分布的伪随机数矩阵
% chi2rnd 可以生成卡方分布的伪随机数矩阵
% trnd 可以生成t分布的伪随机数矩阵
% frnd 可以生成f分布的伪随机数矩阵
% raylrnd   可以生成rayleigh分布的伪随机数矩阵
% lognrnd  对数分布随机抽样

nn=5000; %随机变量个数
mu1=10; std1=1;  %抽样均值和方差
dataa=normrnd(mu1,std1,1,nn); %正态分布抽样
[mua,stda] = normfit(dataa)%正态分布参数估计
fprintf('理论均值10,抽样均值%4.2f 理论方差1,抽样方差%4.2f\n', mua,stda)
%%
%3.连续变量概率密度函数分布和累积函数分布
X=[min(dataa):0.1:max(dataa)];
Y = normpdf(X,mu1,std1);%概率密度函数
p = normcdf(X,mu1,std1);%累积分布函数
% figure(1)
subplot(1,2,1)
plot(X,Y,'r')
xlabel('数据值')
ylabel('对应数据值出现概率')
title('连续变量概率密度函数分布')
% figure(2)
subplot(1,2,2)
plot(X,p,'r')
xlabel('数据值')
ylabel('对应数据值累积概率')
title('连续变量累积概率密度分布')
%%
%4.离散随机变量的概率分布图和直方图
datab=[1 2 4 2 5 3 8 6 5 7 6 3 4 3 4 9 7  3 4 6 5 5 6 4 3 2 5 7  6 5 5 6 6 4 4 4 5 5 1 3 5 ];
mub=mean(datab);
stdb=std(datab);
fprintf('离散变量均值为%4.2f 方程为%4.2f\n', mub,stdb)
[f,xi]=ksdensity(datab);
figure(3)
plot(xi,f,'r')
xlabel('数据值')
ylabel('对应数据值出现概率')
title('随机数据概率密度分布图')
figure(4)
h1 = histogram(datab);
xlabel('数据值')
ylabel('对应数据值出现的次数')
title('随机数据概率直方图')
%%
%5.概率分布模型曲线拟合
%dfittool  %概率分布拟合
save tempresult pd
load  tempresult
rfg=norminv(0.9,pd.mu,pd.sigma);%可靠指标
fprintf('数据小于%4.1f的概率为0.9 \n',rfg)
rfd=norminv(0.1,pd.mu,pd.sigma);%可靠指标
fprintf('数据大于%4.1f的概率为0.9 \n\n',rfd)
p = normcdf([rfg rfd],pd.mu,pd.sigma)%用累积分布函数验证
X1=[min(datab)-2:0.1:max(datab)+2];
Y1 = normpdf(X1,pd.mu,pd.sigma);%概率密度函数
figure(5)
plot(xi,f,'r',X1,Y1,'b.')
xlabel('数据值')
ylabel('对应数据值出现概率')
title('随机数据概率密度分布图')
legend('实际概率分布','拟合概率分布')
%%