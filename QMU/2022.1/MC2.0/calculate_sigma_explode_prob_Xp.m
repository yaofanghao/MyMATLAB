function sigma_explode_prob_Xp = calculate_sigma_explode_prob_Xp(sigma_mean,sigma_variance,prob_p)
%UNTITLED9 此处显示有关此函数的摘要
%   此处显示详细说明
    mu = 0;
    sigma = 1;
    pd = makedist('Normal','mu',mu,'sigma',sigma);
    p = prob_p;
    Up = icdf(pd,p);
    
    sigma_explode_prob_Xp = sqrt(sigma_mean^2 + (Up^2)*(sigma_variance^2));
    
end

