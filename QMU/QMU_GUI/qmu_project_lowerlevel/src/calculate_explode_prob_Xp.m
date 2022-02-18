function [explode_prob_Xp] = calculate_explode_prob_Xp(parameter_miu,prob_p,Std_Dev_sigma)
%UNTITLED5 此处显示有关此函数的摘要
%   此处显示详细说明
    Up = get_Up(prob_p);
    explode_prob_Xp_left = parameter_miu + Up * Std_Dev_sigma;
    explode_prob_Xp_right = parameter_miu - Up * Std_Dev_sigma;
    
    % Change log ： 感度取下限
    explode_prob_Xp = [explode_prob_Xp_left,explode_prob_Xp_right];
end

