function [outputArg1,outputArg2] = plot_result(Xp_upper,Xp_lower,Xpu_upper,Xpu_lower)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
    
x_miritx = [Xp_upper,Xp_lower,Xpu_upper,Xpu_lower,design_safety_X0];
% xlim([0 10])
% ylim([-0.4 0.8])
x = 1:10;
y = x_miritx(1)+zeros(10,1);
plot(x,y,'--b') %蓝色线
hold on;
y = x_miritx(2)+zeros(10,1);
plot(x,y,'--b') %蓝色线
hold on;
y = x_miritx(3)+zeros(10,1);
plot(x,y,'-r') %蓝色线
hold on;
y = x_miritx(4)+zeros(10,1);
plot(x,y,'-r') %蓝色线
hold on;
rectangle('Position',[1 2 5 6])
rymax = 1.1*max(x_miritx);
rymin = 0.9*min(x_miritx);
axis([0 11 rymin rymax]);
end

