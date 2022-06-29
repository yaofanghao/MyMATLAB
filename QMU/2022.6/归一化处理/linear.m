% 线性的归一化方法
clear all;clc;
% 撞击的归一化
miu1 = 1.4667; sigma1 = 0.080838; 
xth1 = 1.1662; xpl1 = 0.71572; x01 = 0.7;
zj_xth = 1; % xth'的值
zj_xpl = xpl1 / zj_xth; % xpl'的值
zj_x0 = x01 / zj_xth; % x0'的值

% 摩擦的归一化
miu2 = 6.24; sigma2 = 0.0627; 
xth2 = 6.0970; xpl2 = 5.6014; x02 = 5.5;
mc_xth = 1;
mc_xpl = xpl2 / mc_xth;
mc_x0 = x02 / mc_xth;

% 多元QMU计算
M = sqrt((zj_xth - zj_x0)^2 + (mc_xth - mc_x0)^2);
U = sqrt((zj_xth - zj_xpl)^2 + (mc_xth - mc_xpl)^2);
Q = M/U;