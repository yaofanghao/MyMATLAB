function [stage_mitrix_i,Stimulus_median_X0,parameter_Step_sized] ...
            = get_stage_mitrix_and_X0(Stimulus_Xi)
%UNTITLED8 此处显示有关此函数的摘要
%   此处显示详细说明
    Stimulus_median_X0 = median( Stimulus_Xi , 'all' );
    temp = Stimulus_Xi - Stimulus_median_X0;
    [~,index] = min(abs(temp));
    Stimulus_median_X0 = Stimulus_Xi(index);
    stage_mitrix_i = zeros();
    
    for i = 1:length(Stimulus_Xi)
        stage_mitrix_i(i) = i - index;
        
    end
    
    % 这里加抛异常
    parameter_Step_sized = abs(Stimulus_Xi(2)-Stimulus_Xi(1));
end

