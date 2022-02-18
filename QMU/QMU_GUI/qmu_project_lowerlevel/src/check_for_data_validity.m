function data_validity = check_for_data_validity(Stimulus_Xi,parameter_M,Distribution_value)
%check_for_data_validity 检查数据有效性
%   此处显示详细说明
    s = size(Stimulus_Xi,1);
    if(Distribution_value)
        if (parameter_M >= 0.25)
            data_validity = 1;
        else
            data_validity = 0;
        end
    else
        if (parameter_M >= 0.3)
            data_validity = 1;
        else
            data_validity = 0;
        end
    end
end

