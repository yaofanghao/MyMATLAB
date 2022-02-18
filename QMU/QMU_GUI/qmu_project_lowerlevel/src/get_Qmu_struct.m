function Qmu_struct = get_Qmu_struct(explode_prob_Xp_QuantileU,explode_prob_Xp,design_safety_X0)
%UNTITLED12 此处显示有关此函数的摘要
%   此处显示详细说明
    Qmu_struct.Qmu_U = abs(max(explode_prob_Xp_QuantileU)-max(explode_prob_Xp));
    Qmu_struct.Qmu_M = abs(design_safety_X0 - max(explode_prob_Xp));
    Qmu_struct.Qmu_Q = Qmu_struct.Qmu_M / Qmu_struct.Qmu_U;
end

