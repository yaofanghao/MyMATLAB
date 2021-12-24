% 需安装AlexNet网络支持包和USB Camera的支持包
clear
camera = webcam; % 连接摄像头
net = alexnet; % 加载网络
while true
  im = snapshot(camera); % 拍一张照片
  image(im); % 展示照片
  im = imresize(im,[227 227]); % resize照片
  label = classify(net,im); % 利用选择的网络分类
  title(char(label)); % 展示分类结果
  drawnow
end

%% 12.23 暂时不知道怎么关闭摄像头


