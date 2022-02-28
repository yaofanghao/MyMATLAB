clc
clear all

xls_dir2 = 'E:\MATLAB\MyMatlab\QMU\MC3.0\moca_QMU\mc_example.xlsx';

mc_data =xlsread(xls_dir2, 'sheet1');

%相关参数设置
zhixindu = 0.95; %置信度
mc_th = 0.65; %安全性设计值

% 计算试验数据的摩擦感度p
explore_num = length(find(mc_data == 1));
mc_prob = explore_num/length(mc_data);

% 估计p的置信区间
% 本质是一元二次方程
mc_Z = get_Up((1-zhixindu) /2);% Za/2
f_a = length(mc_data) +  mc_Z^2;
f_b = -(2*length(mc_data)*mc_prob + mc_Z^2);
f_c = length(mc_data)* (mc_prob^2);
f_delta = sqrt(f_b^2 - 4*f_a*f_c);
mc_p1ow = (-f_b - f_delta) / (2*f_a); %p2
mc_phigh = (-f_b + f_delta) / (2*f_a); %p1

% B类不确定度评定
% 计算sigma_miu
prob_a = (mc_phigh - mc_p1ow) / 2;
mc_v = length(mc_data)-1;
mc_t = tinv(zhixindu, mc_v); %计算t值
mc_ux = prob_a / mc_t;
% 计算sigma_sigma
mc_sigmasigma = mc_ux/sqrt(2* mc_v);
% 计算合成sigma_plow
mc_sigmaplow = sqrt(mc_ux^2 + (mc_t*mc_sigmasigma)^2);

% 计算M
MC_M = mc_p1ow - mc_th;

% 计算U
MC_U = abs(get_Up(zhixindu)) * mc_sigmaplow;

% 计算Q
MC_Q = MC_M / MC_U;





