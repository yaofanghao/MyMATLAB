function [Logistic_row,m,b,n] = calculate_Logistic_row(parameter_M,parameter_A,parameter_n)
%UNTITLED4 此处显示有关此函数的摘要
%   此处显示详细说明
    m = parameter_M;
    b = norm((parameter_A/parameter_n - 1/2),1);
    n = parameter_n;
    
    if (parameter_M >= 0.3)
        Logistic_row = 1.62 * ( parameter_M + 0.029 );
    else
        m = parameter_M;
        b = norm((parameter_A/parameter_n - 1/2),1);
        n = parameter_n;
        prompt = ' M值小于0.3 请查看p(M,b)表输入参数 请查找 ';
        disp(['m = ' m ', b= ' b ', c= ' c]);
        Logistic_row = input(prompt);
    end
end

