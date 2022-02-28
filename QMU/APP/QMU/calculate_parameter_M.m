function parameter_M = calculate_parameter_M(parameter_A,parameter_B,parameter_n)
%UNTITLED3 此处显示有关此函数的摘要
%   此处显示详细说明
    parameter_M = parameter_B / parameter_n - ...
        (parameter_A / parameter_n)^2;
end

