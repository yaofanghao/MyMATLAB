function Up = get_Up(prob_p)
%UNTITLED10 此处显示有关此函数的摘要
%   此处显示详细说明
    mu = 0;
    sigma = 1;
    pd = makedist('Normal','mu',mu,'sigma',sigma);
    p = prob_p;
    Up = icdf(pd,p);
end

