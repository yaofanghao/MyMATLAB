function [parameter_miu, Std_Dev_sigma, explode_prob_Xp, prob, Q] = ...
   shengjiangfa(xint, d, mu, sigma, prob_p, xlim)
%  模拟升降法实验， 输入xint,d,mu,sigma,xlim
%  输出样本miu-hat，sigma-hat，
%  置信水平为prob_p时，爆炸概率响应点Xp的区间，响应概率explode_prob_Xp，Q
%  更新至1.16

n = 50; %样本量
% xci为生成正态随机数
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

% 把模拟数据制作成表格形式
data_example = [xi, v, m, sum_v_m];

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
if data_validity == 0
    fprintf('注意！产生了无效数据\n');
%     error('M小于0.25，生成的数据无效');
end

% 计算mu-hat
parameter_miu = calculate_parameter_miu(...
    Stimulus_median_X0,explode_or_unexplode_data_flag,...
    parameter_A, parameter_n, parameter_Step_sized);

% 计算ρ 
parameter_rou = calculate_parameter_rou(parameter_M,parameter_little_b);
% 计算sigma-hat
Std_Dev_sigma = calculate_Std_Dev_sigma(parameter_rou,parameter_Step_sized);

% 计算置信水平为prob_p时，爆炸概率响应点Xp的区间 
Up = abs(get_Up(prob_p));  %正态的P分位数
[explode_prob_Xp] = calculate_explode_prob_Xp(parameter_miu,prob_p,Std_Dev_sigma);

%响应概率 φ((mu-hat - mu) / sigma)
prob = cdf('Normal',parameter_miu,mu,sigma);
% prob = (length(find(y==1)))/n;

%新定义的Q值计算方法
Q = calculate_parameter_Q(parameter_miu, Up, Std_Dev_sigma, xlim);

end

