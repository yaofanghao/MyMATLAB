%% N次蒙特卡洛试验，调用函数shengjiangfa
% 初始参数设置
num = 100;  % 模拟试验次数N
xint = 1.5;  % 初始刺激量x1
d = 0.1; %步长
mu = 1.4667;  %正态分布参数，（如何设定存疑？）
sigma = 0.0808; 

meann = ones(1, num);  %创建存放均值和标准差的数组
stdn = ones(1, num); 
for j = 1:num
    [a, b] = shengjiangfa(xint, d, mu, sigma);
    meann(:,j) = a;
    stdn(:,j) = b;
end

m = mean(meann);
s = mean(stdn);
fprintf('真实试验mean:%.4f 真实试验std:%.4f \n', mu , sigma)
fprintf('蒙特卡洛模拟mean:%.4f 蒙特卡洛模拟std:%.4f \n', m , s)
disp('success')

% 可视化
xzhou = 1:1:num;
subplot(2,1,1); scatter(xzhou,meann); hold on; plot([0,num],[mu, mu]);
xlabel('试验次数'); ylabel('mean'); axis padded
title('蒙特卡洛模拟升降法');
subplot(2,1,2); scatter(xzhou,stdn); hold on; plot([0,num],[sigma, sigma]);
xlabel('试验次数'); ylabel('std'); axis padded

