%% 2022.1.12- 尝试复现蒙特卡洛升降法
% 
% %% 生成正态随机数
% x = normrnd(1.95, 0.5, 10000, 1)
% % dfittool
% [fp, xp] = ecdf(x)
% ecdfhist(fp, xp, 50)
% hold on
% t = linspace(0, max(x), 100)
% y = normpdf(t, 1.95,0.5)
% plot(t, y, 'g','linewidth',3)
% xlabel('x')
% ylabel('f(x)')
% legend('频率','概率密度')

n = 50; %样本量

% xci为生成正态随机数
% 疑问：mu和sigma的值如何确定？？ 应该来源于真实实验的真值？
xc = normrnd(1.95, 0.5, 1, n);

x = zeros(1,n); %存放xi
x(1) = 2; %设置x1为2
d = 0.5; %设置步长d
y = zeros(1,n); %存放响应值0或1

% 第一次样本试验测试
% if x(1) >= xc(1)
%     y(1) = 1;
%     x(2) = x(1)-d;
% else 
%     y(1) = 0;
%     x(2) = x(1)+d;
% end

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
    fprintf('刺激量个数为%d,数据有效\n', s);
else
    fprintf('刺激量个数为%d,数据无效\n', s);
end

mu1 = mean(x);
sigma1 = std(x);

%响应概率
prob = (length(find(y==1)))/n;

% 可视化
xzhou = 1:1:n;
subplot(2,1,1); scatter(xzhou,xc); hold on; plot(xzhou,x); 
xlabel('样本数'); ylabel('刺激量'); title('临界刺激量和试验刺激量对比'); axis padded

subplot(2,1,2); plot(xzhou,y) 
xlabel('样本数'); ylabel('响应'); title('响应结果'); ylim([0 1]); axis padded

%% 导出数据
xlswrite( 'E:\MATLAB\MyMatlab\QMU\test2022_1\test1.xlsx', xc, 'sheet1', 'A1:AX1')
xlswrite( 'E:\MATLAB\MyMatlab\QMU\test2022_1\test1.xlsx', x,  'sheet1', 'A2:AX2')
xlswrite( 'E:\MATLAB\MyMatlab\QMU\test2022_1\test1.xlsx', y,  'sheet1', 'A3:AX3')
disp('success')