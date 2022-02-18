function [Ni_Stimulus_Xi,explode_or_unexplode_data_flag] = get_Ni_Stimulus_Xi(data_to_process)
%UNTITLED7 此处显示有关此函数的摘要
%   此处显示详细说明
data_to_process = sortrows(data_to_process,1);
sample_N = size(data_to_process,1);
% explode_Ni_Stimulus_Xi 为Xi刺激量下的爆炸次数Ni
explode_Ni_Stimulus_Xi = data_to_process(:,2);

% explode_total_N 为 一次试验中爆炸次数的总和
explode_total_N = sum(explode_Ni_Stimulus_Xi);

% unexplode_Ni_Stimulus_Xi 为 Xi刺激量下的未爆炸次数Ni
unexplode_Ni_Stimulus_Xi = data_to_process(:,3);

% unexplode_total_N 为 一次试验中未爆炸次数的总和
unexplode_total_N = sum(unexplode_Ni_Stimulus_Xi);

if (explode_total_N <= unexplode_total_N)
    Ni_Stimulus_Xi = explode_Ni_Stimulus_Xi;
    explode_or_unexplode_data_flag = -1;
else
    Ni_Stimulus_Xi = unexplode_Ni_Stimulus_Xi;
    explode_or_unexplode_data_flag = 1;
end
end

