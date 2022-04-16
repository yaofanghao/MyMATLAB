clc
clear all

dir = 'image/56.jpg';

img = imread(dir);
imggray = rgb2gray(img);%灰度处理

imggray = imresize(imggray, [256,256]); %图像大小高宽需相等
s = qtdecomp(imggray, 0.2); %四叉树分解，阈值自定
s2 = full(s);

subplot(1,2,1);
imshow(imggray);
subplot(1,2,2);
imshow(s2);