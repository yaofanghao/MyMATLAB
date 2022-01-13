function [meann,stdn] = shengjiangfa(xint, d, mu, sigma)
%SHENGJIANGFA 此处显示有关此函数的摘要
% 模拟升降法实验，输入xint,d,mu,sigma，输出meann,stdn
%   此处显示详细说明

n = 50; %样本量
% xci为生成正态随机数
% 疑问：mu和sigma的值如何确定？？是否应该来源于真实实验的真值？
xc = normrnd(mu, sigma, 1, n);
x = zeros(1,n); %存放xi
x(1) = xint; %初始刺激量
y = zeros(1,n); %存放响应值0或1

for i = 1:n
    if x(i) >= xc(i)
        y(i) = 1;
        x(i+1) = x(i)-d;
    else
        y(i) = 0;
        x(i+1) = x(i)+d;
    end
end
x(:,51) = []; %删除生成的x(51)

% 判定数据有效性：刺激量个数应为4-7个
s = (max(x)-min(x))/d + 1;
if s>=4 || s<=7  
%     fprintf('刺激量个数为%d,数据有效\n', s);
else
    fprintf('刺激量个数为%d,数据无效\n', s); %打印信息检查数据是否有误
end

meann = mean(x);
stdn = std(x);

end

