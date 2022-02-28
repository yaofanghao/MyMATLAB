function parameter_nABMb = calculate_parameter_n_A_B_M_b(Ni_Stimulus_Xi,stage_mitrix_i)
%UNTITLED6 此处显示有关此函数的摘要
%   此处显示详细说明
parameter_nABMb.parameter_n = calculate_parameter_n(Ni_Stimulus_Xi);

parameter_nABMb.parameter_A = calculate_parameter_A(Ni_Stimulus_Xi,stage_mitrix_i);

parameter_nABMb.parameter_B = calculate_parameter_B(Ni_Stimulus_Xi,stage_mitrix_i);

parameter_nABMb.parameter_M = ...
    calculate_parameter_M(parameter_nABMb.parameter_A,parameter_nABMb.parameter_B,parameter_nABMb.parameter_n);

parameter_nABMb.parameter_little_b = ...
    calculate_parameter_little_b(parameter_nABMb.parameter_A,parameter_nABMb.parameter_n);


end

