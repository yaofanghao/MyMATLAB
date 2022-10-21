% 摩擦数据预测
% 机器学习逻辑回归方法
exampleObject = matfile('mc_model.mat');
mcmodel= exampleObject.trainedModel;

mc_table = table(200, 6.35);  %测试数据R=200,压力对数值6.11
% mc_table = table(app.R2Knob.Value, app.NKnob.Value, app.ThetaKnob.Value);
mc_table.Properties.VariableNames{'Var1'} = 'VarName1';
mc_table.Properties.VariableNames{'Var2'} = 'VarName2';
% mc_predict = trainedClassifier.predictFcn(mc_table);

% mc_predict = trainedModel.predictFcn(mc_table); %3.11
mc_predict = mcmodel.predictFcn(mc_table); %4.19

% 神经网络方法
% load ('mcBPmodel.mat');
% p = [200;6.23];
% y = sim(net,p);