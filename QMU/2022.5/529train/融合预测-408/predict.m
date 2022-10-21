clear all
clc

load('BP20.mat')
p = [400;296;66;0];
y = sim(net,p)

%% 
clear all
clc

load('BP20.mat')
p = [500;0;0;1.9];
y = sim(net,p)

