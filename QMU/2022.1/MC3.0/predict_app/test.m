clear all
%% 输入撞击数据的格式

load('zj_model.mat');

zj_table = table(200, 20);  %测试数据R=200,H=20
zj_table.Properties.VariableNames{'Var1'} = 'VarName1';
zj_table.Properties.VariableNames{'Var2'} = 'VarName2';
zj_predict = trainedClassifier.predictFcn(zj_table);

%% 输入摩擦数据的格式
load('mc_net.mat');

%测试数据R=200,P=89,Theta=90
mc_predict = sim(net, [200,  89,  90]'); %注意转置

%% 融合预测概率
alpha = 0.5; %设置撞击概率权重

mix_predict = alpha * zj_predict + (1-alpha) * mc_predict;




