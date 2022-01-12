%--------------------------------------------------------------------------
%  exampA.3-1  基于中位数算法的运动目标检测
%--------------------------------------------------------------------------
% CopyRight：xiezhh

%% 第1种实现方法
%*******************************读取视频数据********************************
% 调用aviread函数读取视频文件
vid = aviread('WalkingMan.avi')
size(vid(1).cdata)
vid(1).colormap

%*****************************数据类型的转换********************************
% 把vid的cdata字段的取值转换成一个240×360×3×80的四维双精度数组IM
IM = cat(4,vid.cdata);
size(IM)
IM = double(IM)/255;


%********************调用median函数求中位数，进行视频图像分割****************
I = median(IM,4);
figure;
imshow(I);                  % 显示背景图像
figure;
imshow(IM(:,:,:,1) - I);    % 显示第1帧中的目标图像



%%  第2种实现方法
WalkManObj = mmreader('WalkingMan.avi');
IM = read(WalkManObj, [1, inf]);
IM = double(IM)/255;
I = median(IM,4);
figure;
imshow(I);                  % 显示背景图像
figure;
imshow(IM(:,:,:,1) - I);    % 显示第1帧中的目标图像


%% 第3种实现方法
WalkManObj = VideoReader('WalkingMan.avi');
IM = WalkManObj.read;
IM = double(IM)/255;
I = median(IM,4);
figure;
imshow(I);                  % 显示背景图像
figure;
imshow(IM(:,:,:,1) - I);    % 显示第1帧中的目标图像