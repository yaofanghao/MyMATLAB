%% 2022.6.29 蒙特卡洛 压力-响应值数据
clc
clear all

% 初始参数设置
R = 400; % 手动选择颗粒直径

num = 50; %样本量

mu = 6.24;  % 压力值取对数 参考常双君专利说明书资料 计算得出
sigma = 0.0627;    

xint = mu;
d = sigma;

% xci为生成正态随机数
xc = normrnd(mu, sigma, 1, num);

x = zeros(1,num); %存放xi
x(1) = xint; %初始刺激量
y = zeros(1,num); %存放响应值0或1

Ri = zeros(1,num); %存放颗粒直径d

for i = 1:num
    Ri(i) = R;
    if x(i) >= xc(i)
        y(i) = 1;
        x(i+1) = x(i)-d;
    else
        y(i) = 0;
        x(i+1) = x(i)+d;
    end
end
x(:,num+1) = []; %删除生成的第num+1个数据

% 存放撞击生成数据集
data_mc = [Ri', x', y']; 

% 可视化
% xzhou = 1:1:num;
% subplot(2,1,1); scatter(xzhou,xc,'+'); hold on; plot(xzhou,x); 
% xlabel('样本数'); ylabel('刺激量'); title('摩擦试验'); axis padded
% 
% subplot(2,1,2); plot(xzhou,y) ;
% xlabel('样本数'); ylabel('响应'); title('响应结果'); ylim([0 1]); axis padded