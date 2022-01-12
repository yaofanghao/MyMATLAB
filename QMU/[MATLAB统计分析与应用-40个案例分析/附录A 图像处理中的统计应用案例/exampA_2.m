%--------------------------------------------------------------------------
%  exampA.2-1  基于K均值聚类的图像分割——灰度图像
%--------------------------------------------------------------------------
% CopyRight：xiezhh

%*******************************读取图像数据********************************
x = imread('coins.png');
whos x
y = double(x(:));


%*************************调用kmeans函数进行聚类分割*************************
startdata = [0; 150];
idpixel = kmeans(y,2,'Start',startdata);
idbw = (idpixel == 2);
result = reshape(idbw, size(x));
imshow(result);



%--------------------------------------------------------------------------
%  exampA.2-2  基于K均值聚类的图像分割——真彩图像
%--------------------------------------------------------------------------

%*******************************读取图像数据********************************
Duck0 = imread('littleduck.jpg');
whos  Duck0
[m,n,k] = size(Duck0);
Duck1 = double(reshape(Duck0, m*n, k));
whos  Duck1


%*************************调用kmeans函数进行聚类分割*************************
startdata = [10 10 200;200 200 10];
idClass = kmeans(Duck1,2,'Start',startdata);
idDuck = (idClass == 1);
result = reshape([idDuck, idDuck, idDuck],[m,n,k]);
Duck2 = Duck0;
Duck2(result) = 0;
figure
imshow(Duck2)