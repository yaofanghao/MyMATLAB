clear
xls_dir = 'E:\MATLAB\MyMatlab\QMU\BP\example.xlsx';
data_example = xlsread(xls_dir, 'sheet1');

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

