function [alert_String,parameter_rou] = ...
    app_request_for_rou_G_H(parameter_M,data_Distribution)
%UNTITLED 计算 rou G H 需要查表进行交互
%   此处显示详细说明
    % 计算rou 对于正态分布和对数正态分布 M > 0.3 公式计算
    % 对于 M <= 0.3 查表
    if(data_Distribution.Distribution_value)
        if (parameter_M > 0.3)
            parameter_rou = 1.62 * ( parameter_M + 0.029 );
            stringM = 'M > 0.3 且为正态分布, rou由公式计算得出,无需查表';
        else
            stringM = '当前分布是正态分布或对数正态分布，请根据上方表格中M,b的值查表 p(M,b), 在右侧输入p(rou)值。';
        end
            stringG = '当前分布是正态分布或对数正态分布，请根据上方表格中p(rou),b的值查表 G(p(rou),b), 在右侧输入G值。';
            stringH = '当前分布是正态分布或对数正态分布，请根据上方表格中p(rou),b的值查表 H(p(rou),b), 在右侧输入H值。';
    else
        stringM = '当前分布是逻辑斯蒂分布或对数逻辑斯蒂分布，请根据上方表格中M,b的值查表 p(M,b), 在右侧输入p(rou)值。';
        stringG = '当前分布是逻辑斯蒂分布或对数逻辑斯蒂分布，请根据上方表格中p(rou),b的值查表 G(p(rou),b), 在右侧输入G值。';
        stringH = '当前分布是逻辑斯蒂分布或对数逻辑斯蒂分布，请根据上方表格中p(rou),b的值查表 H(p(rou),b), 在右侧输入H值。';
    end
    alert_String = {stringM,stringG,stringH};
end

