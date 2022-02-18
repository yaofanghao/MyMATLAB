function output_datasheet = data_display(data_to_process)
%UNTITLED13 此处显示有关此函数的摘要
%   此处显示详细说明
    % 获取刺激矩阵 判定是vi还是mi
    data_to_process = sortrows(data_to_process,1);
    Stimulus_Xi = data_to_process(:,1);
    [Ni_Stimulus_Xi,explode_or_unexplode_data_flag] = ...
        get_Ni_Stimulus_Xi(data_to_process);
    % 获取中位数X0 以及台阶数矩阵,步长d
    [stage_mitrix_i,Stimulus_median_X0,parameter_Step_sized] =...
        get_stage_mitrix_and_X0(Stimulus_Xi);
    % 在末尾添加台阶数
    output_datasheet = [data_to_process,stage_mitrix_i'];
    
end

