function data_validity = check_for_data_validity(parameter_M)
%check_for_data_validity 检查数据有效性
%   此处显示详细说明
    if (parameter_M >= 0.25)
        data_validity = 1;
    else
        data_validity = 0;
%       disp('数据无效');
    end
end

