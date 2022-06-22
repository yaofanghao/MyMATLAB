%% 按下QMU评估按钮 计算QMU
% 撞击Q
xls_dir = 'E:\MATLAB\MyMatlab\QMU\APP\QMU\zj_example.xlsx';
data_example = xlsread(xls_dir, 'sheet1');

[Ni_Stimulus_Xi,explode_or_unexplode_data_flag] = get_Ni_Stimulus_Xi(data_example);

Stimulus_Xi = data_example(:,1); %读取刺激量xi

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

%             %数据有效性判定 1有效 0 无效
%             data_validity = check_for_data_validity(parameter_M);
%             if data_validity == 1
%                 s = '数据有效，请查表输入G,H';
%             else 
%                 s = '数据无效，M<0.25，请重新导入';
%             end
%             
%             % 显示有效性结果
%             app.VaildEditField.Value = s;
%             
%             app.MCNoteEditField.Value = 'G,H,置信水平,安全设计值需自行输入';
%             
%             % 表格数据显示到界面上
%             app.ZJUITable.Data = data_example;

% 计算ρ 
parameter_M = roundn(parameter_M, -2); %保留两位小数
if (parameter_M > 0.3)
parameter_rou = 1.62 * ( parameter_M + 0.029 );
% else
% %                 s = 'M<0.3, 请查表ρ(M,b)，自行输入ρ值！';
% %                 app.MCNoteEditField.Value = s;
%     % 自行给app.rouEditField.Value赋值
%     parameter_rou = app.rouEditField.Value;
end

% 显示计算数据
%             app.MEditField.Value = num2str(parameter_M);
%             app.bEditField.Value = num2str(parameter_little_b);
%             app.rouEditField.Value = num2str(parameter_rou);

% 计算mu-hat
parameter_miu = calculate_parameter_miu(...
Stimulus_median_X0,explode_or_unexplode_data_flag,...
parameter_A, parameter_n, parameter_Step_sized);            

% 计算sigma-hat
Std_Dev_sigma = calculate_Std_Dev_sigma(parameter_rou,parameter_Step_sized);  

% 由自行输入的G H 
% 计算sigma_mean sigma_variance
% parameter_G = app.GEditField.Value;
% parameter_H = app.HEditField.Value;
parameter_H = 1.325;
parameter_G = 1.031;
sigma_mean = parameter_G * Std_Dev_sigma / sqrt(parameter_n);
sigma_variance = parameter_H * Std_Dev_sigma / sqrt(parameter_n);

% 由自行输入的置信水平
% prob = app.UpEditField.Value;
prob = 0.9999;

% 计算sigma_Xp
mu = 0;
sigma = 1;
pd = makedist('Normal','mu',mu,'sigma',sigma);
p = 1 - prob;
Up = icdf(pd,p);
%             p = 0.0001;
sigma_explode_prob_Xp = sqrt(sigma_mean^2 + (Up * sigma_variance)^2 );

% 计算置信区间（上下限）
%             explode_prob_Xp_QuantileU = ...
%                 calculate_explode_prob_Xp_QuantileU(explode_prob_Xp,Confidence_level,sigma_explode_prob_Xp);
%撞击感度取置信下限
% 计算M=Xth-X0
X_th = parameter_miu - abs(Up) * Std_Dev_sigma;
% X_zero = app.X0EditField.Value;
X_zero = 0.7;
ZJ_M = X_th - X_zero;

% 计算U=Xth-Xpl
X_pl = X_th - abs(Up) * sigma_explode_prob_Xp;
ZJ_U = X_th - X_pl; 

% 计算Q
ZJ_Q = ZJ_M / ZJ_U;

%%
% 摩擦Q
% 按下导入摩擦数据按钮
xls_dir2 = 'E:\MATLAB\MyMatlab\QMU\APP\QMU\mc_example.xlsx';

mc_data = xlsread(xls_dir2, 'sheet1');

%相关参数设置
% zhixindu = app.MCalphaEditField.Value;
% mc_th = app.MClimEditField.Value;
zhixindu = 0.95; %置信度
mc_th = 0.956; %安全性设计值

% 计算试验数据的摩擦感度p
explore_num = length(find(mc_data == 1));
mc_prob = explore_num/length(mc_data);

% 估计p的置信区间
% 本质是一元二次方程
mc_Z = get_Up((1-zhixindu) /2);% Za/2
f_a = length(mc_data) +  mc_Z^2;
f_b = -(2*length(mc_data)*mc_prob + mc_Z^2);
f_c = length(mc_data)* (mc_prob^2);
f_delta = sqrt(f_b^2 - 4*f_a*f_c);
mc_p1ow = (-f_b - f_delta) / (2*f_a); %p1
mc_phigh = (-f_b + f_delta) / (2*f_a); %p2
mc_between = [mc_p1ow, mc_phigh];

%             % 显示p和区间 
%             app.P0EditField.Value = num2str(mc_prob);
%             app.MCbetweenEditField.Value = num2str(mc_between);

% B类不确定度评定
% 计算sigma_miu
prob_a = (mc_phigh - mc_p1ow) / 2;
mc_v = length(mc_data)-1;
mc_t = tinv(1-((1-zhixindu) /2), mc_v); %计算t值
mc_ux = prob_a / mc_t;
% 计算sigma_sigma
mc_sigmasigma = mc_ux/sqrt(2* mc_v);
% 计算合成sigma_plow
mc_sigmaplow = sqrt(mc_ux^2 + (mc_t*mc_sigmasigma)^2);

% 计算M
% MC_M = mc_p1ow - mc_th;
MC_M = mc_th - mc_phigh;

% 计算U
% MC_U = abs(get_Up(zhixindu)) * mc_sigmaplow;
MC_U = mc_sigmaplow;

% 计算Q
MC_Q = MC_M ./ MC_U;

% 显示计算数据
% app.MCMEditField.Value = num2str(MC_M);
% app.MCUEditField.Value = num2str(MC_U);
% app.MCQEditField.Value = num2str(MC_Q);
% while MC_Q<=1
%     set(app.MCQLamp, 'Color', [1,0,0])
% end         

%%
% 2.25 测试
% 摩擦Q
clc;clear all;
zhixindu = 0.95; %置信度
mc_th = 1; %安全性设计值
mc_phigh = 0.915;
mc_p1ow = 0.685;

% B类不确定度评定
% 计算sigma_miu
prob_a = (mc_phigh - mc_p1ow) / 2;
mc_v = 50-1;
mc_t = tinv(1-((1-zhixindu) /2), mc_v); %计算t值
mc_ux = prob_a / mc_t;
% 计算sigma_sigma
mc_sigmasigma = mc_ux/sqrt(2* mc_v);
% 计算合成sigma_plow
mc_sigmaplow = sqrt(mc_ux^2 + (mc_t*mc_sigmasigma)^2);

% 计算M
% MC_M = mc_p1ow - mc_th;
MC_M = mc_th - mc_phigh;

% 计算U
% MC_U = abs(get_Up(zhixindu)) * mc_sigmaplow;
MC_U = mc_sigmaplow;

% 计算Q
MC_Q = MC_M ./ MC_U;

% 显示计算数据
% app.MCMEditField.Value = num2str(MC_M);
% app.MCUEditField.Value = num2str(MC_U);
% app.MCQEditField.Value = num2str(MC_Q);
% while MC_Q<=1
%     set(app.MCQLamp, 'Color', [1,0,0])
% end         
