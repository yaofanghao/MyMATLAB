clc
clear all
% 绘制蜘蛛网图
D1 = [1.2, 1.4, 0, 0, 0];
D2 = [1, 1, 0, 0, 0];
P=[D1; D2];

spider_plot(P,...
    'AxesLabels', {'撞击Q', '摩擦Q', '未定义', '未定义', '未定义'},...
    'AxesLimits', [0, 0, 0, 0, 0; 3, 3, 3, 3, 3],... % [min axes limits; max axes limits]
    'AxesPrecision', [1, 1, 1, 1, 1],...
    'AxesInterval', 2,...
    'FillOption', {'on', 'off'},...
    'FillTransparency', [0.2, 0.1] );

legend_str = {'试验感度Q线', '置信阈值线'};
legend(legend_str, 'Location', 'northeast');