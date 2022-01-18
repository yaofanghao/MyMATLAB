clear all

load matlab.mat
%% 预测
predict1 = sim(net, [1.4667 0.0808]')
