function [Ni_Stimulus_Xi,explode_or_unexplode_data_flag] = get_Ni_Stimulus_Xi(data_to_process)
% 比较v和m大小。判定n取vi还是mi
% 输入：要处理的数据（矩阵形式）-------data_to_process
% 输出：n-------Ni_Stimulus_Xi
%       1或-1-------explode_or_unexplode_data_flag
data_to_process = sortrows(data_to_process,1);

explode_Ni_Stimulus_Xi = data_to_process(:,2); %响应vi求和
explode_total_N = sum(explode_Ni_Stimulus_Xi);

unexplode_Ni_Stimulus_Xi = data_to_process(:,3); %不响应mi求和
unexplode_total_N = sum(unexplode_Ni_Stimulus_Xi);

% 比较v和m的大小，取小的作为n
if (explode_total_N <= unexplode_total_N)
    Ni_Stimulus_Xi = explode_Ni_Stimulus_Xi;
    explode_or_unexplode_data_flag = -1;
else
    Ni_Stimulus_Xi = unexplode_Ni_Stimulus_Xi;
    explode_or_unexplode_data_flag = 1;
end
end
