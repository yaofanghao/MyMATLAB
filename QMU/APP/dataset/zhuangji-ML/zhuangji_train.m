clc
clear all

%% 导入数据
xls_dir = 'E:\MATLAB\MyMatlab\QMU\APP\dataset\zhuangji-3.10\test.xlsx';
all = xlsread(xls_dir, 'sheet1');

trainingData= all;
% trainingData = [all(:,1), all(:,2), all(:,3)];

load zj_model_new.mat
% trainingData= zhuangjitrain;

[trainedClassifier, validationAccuracy] = trainClassifier(trainingData);

%% 预测
val = trainingData;
yfit = trainedModel.predictFcn(val);



