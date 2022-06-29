function parameter_n = calculate_parameter_n(Ni_Stimulus_Xi)
%UNTITLED3 此处显示有关此函数的摘要
%   此处显示详细说明
    parameter_n = 0;
    for i = 1:length(Ni_Stimulus_Xi)
        parameter_n = parameter_n + Ni_Stimulus_Xi(i);
    end
end

