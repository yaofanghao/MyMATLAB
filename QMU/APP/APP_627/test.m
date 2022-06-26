%% 按下QMU评估按钮 计算QMU
% 撞击Q
xls_dir = 'E:\MATLAB\MyMatlab\QMU\APP\APP_627\zj_example.xlsx';
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
clear all;clc;
xls_dir2 = 'E:\MATLAB\MyMatlab\QMU\APP\APP_627\mc_example.xlsx';
            
%  输入excel的相对或绝对路径
mc_data = xlsread(xls_dir2, 'sheet1');

 % 比较v和m大小。判定n取vi还是mi
[Ni_Stimulus_Xi2,~] = get_Ni_Stimulus_Xi(mc_data);

Stimulus_Xi2 = mc_data(:,1); %读取刺激量xi

% 获取中位数X0 以及台阶数矩阵,步长d
% [~,~,~] =...
%     get_stage_mitrix_and_X0(Stimulus_Xi2);
stage_mitrix_i2 = [0,1,2,3,4]; % 因为摩擦的 台阶数i 的计算方法和撞击不同，这里强制修改

% 计算A B M b 的模块
nABMb2 = calculate_parameter_n_A_B_M_b(Ni_Stimulus_Xi2,stage_mitrix_i2);
parameter_n2 = nABMb2.parameter_n;
parameter_A2 = nABMb2.parameter_A;
parameter_B2 = nABMb2.parameter_B;
parameter_M2 = nABMb2.parameter_M;
parameter_little_b2 = nABMb2.parameter_little_b;      

%             %数据有效性判定 1有效 0 无效
%             data_validity = check_for_data_validity(parameter_M);
%             if data_validity == 1
%                 s = '数据有效，请查表输入G,H';
%             else 
%                 s = '数据无效，M<0.25，请重新导入';
%             end
%             % 显示有效性结果
%             app.VaildEditField.Value = s;

% 表格数据显示到界面上
% app.MCUITable.Data = app.mc_data;

%---------------------------------------------------------------------            

% 计算ρ 
parameter_M2 = roundn(parameter_M2, -2); %保留两位小数
if (parameter_M2 > 0.3)
    parameter_rou2 = 1.62 * ( parameter_M2 + 0.029 );
else
%                 s = 'M<0.3, 请查表ρ(M,b)，自行输入ρ值！';
%                 app.NoteEditField.Value = s;
    % 自行给app.rouEditField.Value赋值
    parameter_rou2 = rou2EditField.Value;
end

% 按下导入摩擦数据按钮
[~,explode_or_unexplode_data_flag2] = get_Ni_Stimulus_Xi(mc_data);

% 获取中位数X0 以及台阶数矩阵,步长d
[~,Stimulus_median_X02,parameter_Step_sized2] =...
    get_stage_mitrix_and_X0(Stimulus_Xi2);

% 计算mu-hat
parameter_miu2 = calculate_parameter_miu(...
    Stimulus_median_X02,explode_or_unexplode_data_flag2,...
    parameter_A2, parameter_n2, parameter_Step_sized2);            

% 计算sigma-hat
Std_Dev_sigma2 = calculate_Std_Dev_sigma(parameter_rou2,parameter_Step_sized2);  

% 由自行输入的G H 
% 计算sigma_mean sigma_variance
parameter_G2 = 0.961;
parameter_H2 = 1.596;
sigma_mean2 = parameter_G2 * Std_Dev_sigma2 / sqrt(parameter_n2);
sigma_variance2 = parameter_H2 * Std_Dev_sigma2 / sqrt(parameter_n2);

% 由自行输入的置信水平
% prob2 = str2double(MCalphaEditField.Value);
prob2 = 0.9999;
% 计算sigma_Xp
mu2 = 0;
sigma2 = 1;
pd2 = makedist('Normal','mu',mu2,'sigma',sigma2);
p2 = 1 - prob2;
Up2 = icdf(pd2,p2);
%             p = 0.0001;
sigma_explode_prob_Xp2 = sqrt(sigma_mean2^2 + (Up2 * sigma_variance2)^2 );

% 计算置信区间（上下限）
% 取置信下限
% 计算M=Xth-X0
X_th2 = parameter_miu2 - abs(Up2) * Std_Dev_sigma2;

%             X_zero = str2double(app.X0EditField.Value);
% 显示安全性设计值X0，若满足不等式大于1，X0小于X_pl即可

X_pl2 = X_th2 - abs(Up2) * sigma_explode_prob_Xp2;             

% if app.DropDown.Value == char('自动计算')
    mc_th = X_pl2 - 0.1;
% else
%     mc_th = str2double(app.MClimEditField.Value); % 可以自定义修改安全性设计值
% end

MC_M = X_th2 - mc_th;

% 计算U
MC_U = X_th2 - X_pl2;


% 计算Q
MC_Q = MC_M ./ MC_U;
            