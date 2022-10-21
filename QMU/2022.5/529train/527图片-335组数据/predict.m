clear all
%% classification learner
% 输入数据的格式-代码未成功 5.27

exampleObject = matfile('mix_model.mat');
mixmodel= exampleObject.trainedModel;

mix_table = table(400, 470, 90, 1.2);  %测试数据R=200,H=20
mix_table.Properties.VariableNames{'Var1'} = 'MPa';
mix_table.Properties.VariableNames{'Var2'} = 'VarName1';
mix_table.Properties.VariableNames{'Var3'} = 'VarName3';
mix_table.Properties.VariableNames{'Var4'} = 'VarName4';
mix_predict = mixmodel.predictFcn(mix_table);
