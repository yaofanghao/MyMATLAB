clc;clear;
data_example = xlsread('./datapaper.xlsx');
Stimulus_Xi = data_example(:,1);
data_Distribution = Dataset_Property_Class(1);
% 获取刺激矩阵 判定是vi还是mi
[Ni_Stimulus_Xi,explode_or_unexplode_data_flag] = ...
    get_Ni_Stimulus_Xi(data_example);
% 获取中位数X0 以及台阶数矩阵,步长d
[stage_mitrix_i,Stimulus_median_X0,parameter_Step_sized] =...
    get_stage_mitrix_and_X0(Stimulus_Xi);
% test data
stage_mitrix_i = [0,1,2,3];
Stimulus_median_X0 = 23;

nABMb = calculate_parameter_n_A_B_M_b(Ni_Stimulus_Xi,stage_mitrix_i);

parameter_n = nABMb.parameter_n;

parameter_A = nABMb.parameter_A;

parameter_B = nABMb.parameter_B;

parameter_M = nABMb.parameter_M;

parameter_little_b = nABMb.parameter_little_b;
% 数据有效性判定
data_validity = check_for_data_validity(Stimulus_Xi,parameter_M,data_Distribution);

parameter_miu = calculate_parameter_miu(...
    Stimulus_median_X0,explode_or_unexplode_data_flag,...
    parameter_A, parameter_n, parameter_Step_sized);
    
[parameter_rou,parameter_G,parameter_H] = ...
    request_for_the_rou_and_G_H(parameter_M,parameter_little_b,data_Distribution);
% parameter_rou = 0.801 G = 1.622 H = 1.643
Std_Dev_sigma = calculate_Std_Dev_sigma(parameter_rou,parameter_Step_sized);

[sigma_mean,sigma_variance] = ...
    calculate_sigma_mean_variance(parameter_G,parameter_H,parameter_n,Std_Dev_sigma);
% 先取0.9 用于测试
prob_p = 1e-6;   Up = abs(get_Up(prob_p));
[explode_prob_Xp] = calculate_explode_prob_Xp(parameter_miu,prob_p,Std_Dev_sigma);

sigma_explode_prob_Xp = ...
    calculate_sigma_explode_prob_Xp(sigma_mean,sigma_variance,prob_p);
Confidence_level = 0.9999;
explode_prob_Xp_QuantileU = ...
    calculate_explode_prob_Xp_QuantileU(explode_prob_Xp(2),Confidence_level,sigma_explode_prob_Xp);
design_safety_X0 = 26;
Qmu = get_Qmu_struct(explode_prob_Xp_QuantileU,explode_prob_Xp,design_safety_X0);




