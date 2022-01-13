%% N次蒙特卡洛试验，调用函数shengjiangfa
% 初始参数设置
num = 100;  % 模拟试验次数N
xint = 2;  % 初始刺激量x1
d = 0.5; %步长
mu = 1.95; %正态分布参数，（如何设定存疑？）
sigma = 0.5; 

meann = ones(1, num);  %创建存放均值和标准差的数组
stdn = ones(1, num); 
for j = 1:num
    [a, b] = shengjiangfa(xint, d, mu, sigma);
    meann(:,j) = a;
    stdn(:,j) = b;
end

mean(meann)
mean(stdn)

disp('success')

% 可视化
xzhou = 1:1:num;
subplot(2,1,1); scatter(xzhou,meann); hold on; plot([0,num],[1.95,1.95]);
xlabel('试验次数'); ylabel('mean'); axis padded
title('蒙特卡洛模拟升降法');
subplot(2,1,2); scatter(xzhou,stdn); hold on; plot([0,num],[0.5,0.5]);
xlabel('试验次数'); ylabel('std'); axis padded

