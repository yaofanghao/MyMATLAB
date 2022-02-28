clear all
%% 输入撞击数据的格式
load zj_model.mat;

zj_table = table(200, 20);  %测试数据R=200,H=20
zj_table.Properties.VariableNames{'Var1'} = 'VarName1';
zj_table.Properties.VariableNames{'Var2'} = 'VarName2';
zj_predict = trainedModel1.predictFcn(zj_table);

%% 输入摩擦数据的格式
clear all
load mc_model.mat;

mc_table = table(200, 270, 30);  %测试数据R=200,压力270,摆角30
mc_table.Properties.VariableNames{'Var1'} = 'VarName1';
mc_table.Properties.VariableNames{'Var2'} = 'VarName2';
mc_table.Properties.VariableNames{'Var3'} = 'VarName3';
% mc_predict = trainedClassifier.predictFcn(mc_table);
mc_predict = trainedModel.predictFcn(mc_table);

%% 融合预测概率
alpha = 0.5; %设置撞击概率权重

mix_predict = alpha * zj_predict + (1-alpha) * mc_predict;




