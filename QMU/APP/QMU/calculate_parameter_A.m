function parameter_A = calculate_parameter_A(Ni_Stimulus_Xi,stage_mitrix_i)
%UNTITLED3 此处显示有关此函数的摘要
%   此处显示详细说明
    parameter_A = 0;
    for i = 1:length(Ni_Stimulus_Xi)
        parameter_A = parameter_A + stage_mitrix_i(i) * Ni_Stimulus_Xi(i);
    end
end

