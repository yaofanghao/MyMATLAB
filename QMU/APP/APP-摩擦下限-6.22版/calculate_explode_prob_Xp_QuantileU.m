function explode_prob_Xp_QuantileU = ...
    calculate_explode_prob_Xp_QuantileU(explode_prob_Xp,Confidence_level,sigma_explode_prob_Xp)
%UNTITLED11 此处显示有关此函数的摘要
%   此处显示详细说明
    Up = get_Up(Confidence_level);
    explode_prob_Xp_QuantileU_left = explode_prob_Xp - Up * sigma_explode_prob_Xp;
    explode_prob_Xp_QuantileU_right = explode_prob_Xp + Up * sigma_explode_prob_Xp;
    explode_prob_Xp_QuantileU = [explode_prob_Xp_QuantileU_left,explode_prob_Xp_QuantileU_right];
end

