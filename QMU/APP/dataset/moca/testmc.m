%% 输入摩擦数据的格式
load mc_model.mat;

mc_table = table(200, 269, 67);  %测试数据R=200,压力270,摆角30
mc_table.Properties.VariableNames{'Var1'} = 'VarName1';
mc_table.Properties.VariableNames{'Var2'} = 'VarName2';
mc_table.Properties.VariableNames{'Var3'} = 'VarName3';
% mc_predict = trainedClassifier.predictFcn(mc_table);
mc_predict = trainedModel.predictFcn(mc_table);
