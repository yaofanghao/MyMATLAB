function parameter_B = calculate_parameter_B(Ni_Stimulus_Xi,stage_mitrix_i)
%UNTITLED3 此处显示有关此函数的摘要
%   此处显示详细说明
    parameter_B = 0;
    for i = 1:length(Ni_Stimulus_Xi)
        parameter_B = parameter_B + stage_mitrix_i(i) * ...
                      stage_mitrix_i(i) * Ni_Stimulus_Xi(i);
    end
end