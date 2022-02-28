%% 输入撞击数据的格式
load zj_model.mat;

zj_table = table(200, 20);  %测试数据R=200,H=20
zj_table.Properties.VariableNames{'Var1'} = 'VarName1';
zj_table.Properties.VariableNames{'Var2'} = 'VarName2';
zj_predict = trainedModel1.predictFcn(zj_table);
