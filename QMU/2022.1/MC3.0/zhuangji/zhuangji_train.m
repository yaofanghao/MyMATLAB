clc
clear all
%导入数据
% xls_dir = 'C:\Users\姚方浩\Desktop\MC3.0\zhuangji\zhuangji.xlsx';
% all = xlsread(xls_dir, 'sheet1');
% 
% trainingData= all;
% % trainingData = [all(:,1), all(:,2), all(:,3)];

load model.mat
trainingData= zhuangjitrain;

[trainedClassifier, validationAccuracy] = trainClassifier(trainingData);

%% 预测
val = zhuangjivalid;
yfit = trainedClassifier.predictFcn(val);



