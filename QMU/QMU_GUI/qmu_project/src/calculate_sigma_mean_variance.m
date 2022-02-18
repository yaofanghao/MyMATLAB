function [sigma_mean,sigma_variance] = ...
    calculate_sigma_mean_variance(parameter_G,parameter_H,parameter_n,Std_Dev_sigma)
%UNTITLED8 此处显示有关此函数的摘要
%   此处显示详细说明
    sigma_mean = parameter_G * Std_Dev_sigma / sqrt(parameter_n);
    sigma_variance = parameter_H * Std_Dev_sigma / sqrt(parameter_n);
end

