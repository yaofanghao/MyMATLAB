function result_table = generate_result_table(Data_Whole_Structure)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
result_table = cell(30,1);
result_table{1,1} = Data_Whole_Structure.Distribution;
result_table{2,1} = Data_Whole_Structure.Distribution_value;
result_table{3,1} = Data_Whole_Structure.Qmu_Q;
result_table{4,1} = Data_Whole_Structure.Qmu_M;
result_table{5,1} = Data_Whole_Structure.Qmu_U;
result_table{6,1} = Data_Whole_Structure.explode_or_unexplode_data_flag;
result_table{7,1} = Data_Whole_Structure.Stimulus_median_X0;
result_table{8,1} = Data_Whole_Structure.parameter_Step_sized;
result_table{9,1} = Data_Whole_Structure.parameter_n;
result_table{10,1} = Data_Whole_Structure.parameter_A;
result_table{11,1} = Data_Whole_Structure.parameter_B;
result_table{12,1} = Data_Whole_Structure.parameter_M;
result_table{13,1} = Data_Whole_Structure.parameter_little_b;
result_table{14,1} = Data_Whole_Structure.data_validity;
result_table{15,1} = Data_Whole_Structure.prob_p;
result_table{16,1} = Data_Whole_Structure.Up;
result_table{17,1} = Data_Whole_Structure.parameter_miu;
result_table{18,1} = Data_Whole_Structure.parameter_rou;
result_table{19,1} = Data_Whole_Structure.parameter_G;
result_table{20,1} = Data_Whole_Structure.parameter_H;
result_table{21,1} = Data_Whole_Structure.Std_Dev_sigma;
result_table{22,1} = Data_Whole_Structure.sigma_mean;
result_table{23,1} = Data_Whole_Structure.sigma_variance;
result_table{24,1} = Data_Whole_Structure.Confidence_level;
result_table{25,1} = Data_Whole_Structure.design_safety_X0;
result_table{26,1} = Data_Whole_Structure.sigma_explode_prob_Xp;

result_table{27,1} = Data_Whole_Structure.explode_prob_Xp_upper;
result_table{28,1} = Data_Whole_Structure.explode_prob_Xp_lower;

result_table{29,1} = Data_Whole_Structure.explode_prob_Xp_QuantileU_upper;
result_table{30,1} = Data_Whole_Structure.explode_prob_Xp_QuantileU_lower;

end

