function parameter_little_b = calculate_parameter_little_b(parameter_A,parameter_n)
%calculate_parameter_little_b 计算小b， 取小数部分四舍五入一位
%   此处显示详细说明
    temp = abs(parameter_A / parameter_n - 0.5);
    temp_mis = temp - floor(temp);
    temp_mis = round(temp_mis,1);
    if temp_mis <= 0.5
        parameter_little_b = temp_mis;
    else
        parameter_little_b = 1 - temp_mis;
    end
end

