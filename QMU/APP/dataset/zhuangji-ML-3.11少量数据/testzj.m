%% 输入撞击数据的格式
load zj_model_new.mat;

zj_table = table(400, 1.4);  %测试数据R,H
zj_table.Properties.VariableNames{'Var1'} = 'VarName1';
zj_table.Properties.VariableNames{'Var2'} = 'VarName2';
zj_predict = trainedModel.predictFcn(zj_table);

%% 批量预测撞击
% 未能成功运行 3.10
load zj_model_new.mat;
% 先导入test.xlsx

num = 375;
zj_predict = zeros(1,num);
for i = 1:num
    
zj_table = table(test(i,1), test(i,2));  % R,H=
zj_table.Properties.VariableNames{'Var1'} = 'VarName1';
zj_table.Properties.VariableNames{'Var2'} = 'VarName2';
zj_predict(i) = trainedModel.predictFcn(zj_table);

end


