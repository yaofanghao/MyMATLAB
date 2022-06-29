clc
clear all

%% 2022.1.12 尝试蒙特卡洛升降法--单次
% 更新至1.18

% 初始参数设置
n = 50; %样本量
xint = 1.5;  % 初始刺激量x1
d = 0.1; %步长
mu = 1.4667;  %？正态分布参数（依据真实试验作为先验知识）
sigma = 0.0808; 
prob_p = 1e-4;  %置信度水平
xlim = 0.7; %？根据先验知识取临界安全值（由特性落高安全性设计值公式推出） 

% xci为生成正态随机数
xc = normrnd(mu, sigma, 1, n);
x = zeros(1,n); %存放xi
x(1) = xint; %初始刺激量
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
    fprintf('刺激量个数为%s,数据有效\n', s);
else
    error('刺激量个数为%s,数据无效\n', s);
end

%% 导出数据到test1.xlsx （非必需）
xlswrite_dir = 'E:\MATLAB\MyMatlab\QMU\MC2.0\test1.xlsx';
xlswrite( xlswrite_dir, xc, 'sheet1', 'A1:AX1')
xlswrite( xlswrite_dir, x,  'sheet1', 'A2:AX2')
xlswrite( xlswrite_dir, y,  'sheet1', 'A3:AX3')
disp('success')
% xlswrite('filename',' ','Sheet1')

%% 转化为实验数据型的表格
% 参考格式：example.xlsx
table = tabulate(x); %读取xi的数值和频数
xi = table(:,1); % 刺激量xi
sum_v_m = table(:,2); % xi的频数

%创建响应和不响应的列
v = zeros(length(xi),1);
m = zeros(length(xi),1); 

%统计刺激量为xi(j)时，响应和不响应个数
for j = 1:length(xi)
    for i = 1:n
        if  x(i) == xi(j)  
            if y(i)==1
                v(j) = v(j) + 1;
            else 
                m(j) = m(j) + 1;
            end
        end
    end
end

%% 上一模块数据处理后导出到test2.xlsx （非必需）
xls_dir = 'E:\MATLAB\MyMatlab\QMU\MC2.0\test2.xlsx';
xlswrite( xls_dir, xi, 'sheet1', 'A1')
xlswrite( xls_dir, v,  'sheet1', 'B1')
xlswrite( xls_dir, m,  'sheet1', 'C1')
xlswrite( xls_dir, sum_v_m,  'sheet1', 'D1')
disp('success')

% 读取表格信息
% data_example = xlsread( xls_dir, 'sheet1');
%% 计算爆炸临界阈值mu-hat和方差sigma-hat

% 把模拟数据制作成表格形式
data_example = [xi, v, m, sum_v_m];
data_Distribution = 1; %正态或对数正态分布
Stimulus_Xi = data_example(:,1); %读取刺激量xi

% 比较v和m大小。判定n取vi还是mi
[Ni_Stimulus_Xi,explode_or_unexplode_data_flag] = get_Ni_Stimulus_Xi(data_example);

% 获取中位数X0 以及台阶数矩阵,步长d
[stage_mitrix_i,Stimulus_median_X0,parameter_Step_sized] =...
    get_stage_mitrix_and_X0(Stimulus_Xi);

% 计算A B M b 的模块
nABMb = calculate_parameter_n_A_B_M_b(Ni_Stimulus_Xi,stage_mitrix_i);
parameter_n = nABMb.parameter_n;
parameter_A = nABMb.parameter_A;
parameter_B = nABMb.parameter_B;
parameter_M = nABMb.parameter_M;
parameter_little_b = nABMb.parameter_little_b;

% 数据有效性判定 1有效 0 无效
data_validity = check_for_data_validity(parameter_M);

% 计算mu-hat
parameter_miu = calculate_parameter_miu(...
    Stimulus_median_X0,explode_or_unexplode_data_flag,...
    parameter_A, parameter_n, parameter_Step_sized);

% 计算ρ 
parameter_rou = calculate_parameter_rou(parameter_M,parameter_little_b);
if parameter_rou==0
   error('M小于0.25，数据无效！') 
end

% 计算sigma-hat
Std_Dev_sigma = calculate_Std_Dev_sigma(parameter_rou,parameter_Step_sized);

%% 估计给定置信度下的响应概率（Xp）区间和新定义Q值计算方法，作为模型训练集用的输出
% 计算置信水平为prob_p时，爆炸概率响应点Xp的区间
prob_p = 1e-4;   
Up = abs(get_Up(prob_p));  %正态的P分位数
[explode_prob_Xp] = calculate_explode_prob_Xp(parameter_miu,prob_p,Std_Dev_sigma);

%响应概率 φ((mu-hat - mu) / sigma)
prob = cdf('Normal',parameter_miu,mu,sigma);
% prob = (length(find(y==1)))/n;

%新定义的Q值计算方法
Q = calculate_parameter_Q(parameter_miu, Up, Std_Dev_sigma, xlim);

%% 模拟试验结果打印
xzhou = 1:1:n;
subplot(2,2,1); scatter(xzhou,xc,'+'); hold on; plot(xzhou,x); 
xlabel('样本数'); ylabel('刺激量'); title('临界刺激量和试验刺激量对比'); axis padded

subplot(2,2,2); plot(xzhou,y) ;
xlabel('样本数'); ylabel('响应'); title('响应结果'); ylim([0 1]); axis padded

subplot(2,2,3); plot([0, n],[prob, prob],'blue'); hold on; plot([0, n],[0.5, 0.5],'red');
xlabel('样本数'); ylabel('Xp'); ylim([0 1]); legend('试验的Xp','0.5概率线'); 
% axis padded
title('Xp'); 

subplot(2,2,4); plot([0, n],[Q,Q], 'blue'); hold on; plot([0, n],[1,1],'red')
xlabel('样本数'); ylabel('Q'); ylim([-1 4]); legend('试验的Q值','Q=1边界线'); 
% axis padded
title('Q值'); 

fprintf('蒙特卡洛模拟mu-hat:%.4f sigma-hat:%.4f \n', parameter_miu, Std_Dev_sigma)
fprintf('响应概率：%.2f \n估计Q值：%.4f \n置信区间:', prob, Q)
disp(explode_prob_Xp)
