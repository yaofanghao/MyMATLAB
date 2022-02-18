function [parameter_rou,parameter_G,parameter_H] = ...
    request_for_the_rou_and_G_H(parameter_M,parameter_little_b,data_Distribution)
%UNTITLED 计算 rou G H 需要查表进行交互
%   此处显示详细说明
    % 计算rou 对于正态分布和对数正态分布 M > 0.3 公式计算
    % 对于 M <= 0.3 查表
    if(data_Distribution.Distribution_value)
        if (parameter_M > 0.3)
            parameter_rou = 1.62 * ( parameter_M + 0.029 );
        else
            string_cur_dis = ['当前分布是 ', data_Distribution.Distribution];
            disp(string_cur_dis);
            string = '请查表 p(M,b), M = %6.3f b = %6.3f \n';
            fprintf(string,parameter_M,parameter_little_b);
            parameter_rou = input('请输入 p(M,b) 值');
        end
            string_cur_dis = ['当前分布是 ', data_Distribution.Distribution];
            disp(string_cur_dis);
            string = '请查表 G(p(rou),b), p(rou) = %6.3f b = %6.3f \n';
            fprintf(string,parameter_rou,parameter_little_b);
            parameter_G = input('请输入 G(p(rou),b) 值:');

            string = '请查表 H(p(rou),b), p(rou) = %6.3f b = %6.3f \n';
            fprintf(string,parameter_rou,parameter_little_b);
            parameter_H = input('请输入 H(p(rou),b) 值:');
    else
        % 计算rou 对于逻辑斯蒂分布和对数逻辑斯蒂分布
        string_cur_dis = ['当前分布是 ', data_Distribution.Distribution];
        disp(string_cur_dis);
        string = '请查表 p(M,b), M = %6.3f b = %6.3f \n';
        fprintf(string,parameter_M,parameter_little_b);
        parameter_rou = input('请输入 p(M,b) 值');
        string = '请查表 G(p(rou),b), p(rou) = %6.3f b = %6.3f \n';
        fprintf(string,parameter_rou,parameter_little_b);
        parameter_G = input('请输入 G(p(rou),b) 值');
        
        string = '请查表 H(p(rou),b), p(rou) = %6.3f b = %6.3f \n';
        fprintf(string,parameter_rou,parameter_little_b);
        parameter_H = input('请输入 H(p(rou),b) 值');
    end
end

