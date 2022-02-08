close all;
clear all;
clc; 

%% MATLAB自带函数matchFeatures的使用示例
% indexPairs = matchFeatures(features1,features2); 
% [indexPairs,matchmetric] = matchFeatures(features1,feature2); 
% [indexPairs,matchmetric]=matchFeatures(features1,feature2,Name,Value); 
% features1,features2是提取出的特征描述子
% indexPairs为nx2的向量，即互相匹配的n对点所对应的坐标索引
% matchmetric为匹配后的特征描述子之间的测度值
% Name为用一对单引号包含的字符串，Value为对应Name的值。

I1 = rgb2gray(imread('left.png')); 
I2 = rgb2gray(imread('right.png'));

points1 = detectHarrisFeatures(I1); 
points2 = detectHarrisFeatures(I2); 

[features1,valid_points1] = extractFeatures(I1,points1); 
[features2,valid_points2] = extractFeatures(I2,points2); 

indexPairs = matchFeatures(features1,features2); 

matchedPoints1 = valid_points1(indexPairs(:,1)); 
matchedPoints2 = valid_points2(indexPairs(:,2)); 

figure; 
showMatchedFeatures(I1,I2,matchedPoints1,matchedPoints2);
