clc
clear all
% 2022.1.12 尝试蒙特卡洛升降法--N次 调用函数shengjiangfa
% 更新至1.16

% 初始参数设置
num = 20;  % 模拟试验次数N
xint = 1.5;  % 初始刺激量x1
d = 0.1; %步长
mu = 1.4667;  %？正态分布参数（依据真实试验作为先验知识）
sigma = 0.0808; 
prob_p = 1e-4;  %置信度水平，这里取0.9999
x_lim = 0.7; %？根据先验知识取临界安全值（由特性落高安全性设计值公式推出） 

mu_hatn = zeros(num, 1 );  %创建存放N次试验的均值 标准差 Xp Q值
sigma_hatn = zeros(num, 1); 
prob_n  = zeros(num, 1);
Q_hatn = zeros(num,1);

for j = 1:num
    [mu_hat, sigma_hat, Xp, prob, Q_hat] = ...
   shengjiangfa(xint, d, mu, sigma, prob_p, x_lim); %shengjiangfa函数3.0【改进版】
    
    mu_hatn(j,:) = mu_hat;
    sigma_hatn(j,:) = sigma_hat;
    prob_n(j,:) = prob;
    Q_hatn(j,:) = Q_hat;

    fprintf('mu-hat:%.4f sigma-hat:%.4f ', mu_hat,sigma_hat)
    fprintf('prob：%.4f Q_hat：%.4f Xp:', prob, Q_hat)
    disp(Xp)
end

fail_num = 0;
less_than_one = 0;
for j = 1:num
     if Q_hatn(j,:)==Inf
        fprintf('第%d组的数据无效,M小于0.25\n',j)
        fail_num = fail_num +1;
    elseif Q_hatn(j,:)<1 
        fprintf('第%d组的数据Q值小于1 \n',j)
        less_than_one = less_than_one +1;         
    end
end

fprintf('-----------------------\n')
fprintf('本模拟试验共%d组数据,其中%d组数据有效。\n',num, num-fail_num)
fprintf('%d组数据满足置信水平%.4f时，',num - fail_num - less_than_one, 1-prob_p)
fprintf('爆炸概率不高于%.4f的安全性要求。\n', prob_p)
fprintf('-----------------------\n')

% 制作数据集
make_dataset = [mu_hatn, sigma_hatn, prob_n, Q_hatn];
%% 可视化
xzhou = 1:1:num;
subplot(2,2,1); scatter(xzhou,mu_hatn,'blue'); hold on; plot([0,num],[mu, mu]);
ylabel('mu-hat'); axis padded
title('mu-hat');
subplot(2,2,2); scatter(xzhou,sigma_hatn,'blue'); hold on; plot([0,num],[sigma, sigma]);
ylabel('sigma-hat'); axis padded
title('sigma-hat');
subplot(2,2,3); scatter(xzhou,prob_n,'red'); hold on; plot([0, num],[0.5, 0.5]);
ylabel('Xp'); ylim([0 1])
title('Xp'); 
subplot(2,2,4); scatter(xzhou,Q_hatn,'red'); hold on; plot([0, num],[1,1]);
xlabel('试验次数'); ylim([-1 4]);
%legend('试验的Q值','边界Q=1');
title('Q');


