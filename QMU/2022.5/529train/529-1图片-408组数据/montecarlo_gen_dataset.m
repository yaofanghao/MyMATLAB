
%% 2022.5.27 尝试蒙特卡洛生成融合数据
clc
clear all

% 初始参数设置
% R = 200; % 手动选择颗粒直径
% R = 300;
R = 400;
% R = 500;

num = 20; %样本量

if R == 200 
    mu = 82.5;  % 临界压力的正态分布参数（？依据仿真试验产生数据作为先验知识）
    sigma = 7.5; 
    
end
if R == 300
    mu = 90;
    sigma = 5;
end
if R == 400
    mu =0.3;
    sigma = 0.15;
end
if R == 500
   mu = 96;
   sigma = 1;
end
xint = mu;  % 初始刺激量x1
d = sigma;  %步长

% xci为生成正态随机数
xc = normrnd(mu, sigma, 1, num);

x = zeros(1,num); %存放xi
x(1) = xint; %初始刺激量
y = zeros(1,num); %存放响应值0或1

Ri = zeros(1,num); %存放颗粒直径d
thetai = zeros(1,num);

for i = 1:num
    Ri(i) = R;
    thetai(i) = 90; %摆角
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
data_mc = [Ri', x', thetai', y']; 

% % 可视化
% xzhou = 1:1:num;
% subplot(2,1,1); scatter(xzhou,xc,'+'); hold on; plot(xzhou,x); 
% xlabel('样本数'); ylabel('刺激量'); title('摩擦试验'); axis padded
% 
% subplot(2,1,2); plot(xzhou,y) ;
% xlabel('样本数'); ylabel('响应'); title('响应结果'); ylim([0 1]); axis padded