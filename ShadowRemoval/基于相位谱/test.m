%% Murali 和 Govindan - Removal of Shadows from a Single Image.pdf
clc
clear all

dir = 'image/56.jpg';
img = imread(dir);

% 1. 将RGB图像转换为LAB颜色空间。
% L = mean2(img); % 平均亮度
R = img(:,:,1);
G = img(:,:,2);
B = img(:,:,3);

Acolor = R./G;
Bcolor = R./B;

% hsv空间
[m,n,k] = size(img); %读取图片大小
hsv = rgb2hsv(img); %颜色空间转换
H = hsv(:,:,1); % 色调
S = hsv(:,:,2); % 饱和度
V = hsv(:,:,3); % 亮度

AverB = mean2(Bcolor);
AverV = mean2(V);
%%  2.将L和B平面上数值较低的像素分类为阴影像素，其他像素为非阴影像素（基于一个阈值）。

for i = 1:m %遍历每一个像素点，可以根据需要选择自己需要处理的区域
    for j = 1: n
        if (hsv(i,j,3) < AverV) && (img(i,j,3) < AverB)
            hsv(i,j,3) = AverV*2;
%             img(i,j,3) = AverB;
        end
        
    end
end

rgb1 = hsv2rgb(hsv); %转为RGB，进行显示

subplot(1,2,1);imshow(img)
subplot(1,2,2);imshow(rgb1)

%% 3.分别定位阴影区域。

%% 4. 对于每个阴影区域，做以下工作 
% 4.1. 计算阴影区域内的R、G和B通道的平均数值。
% 4.2. 计算该区域外的非阴影区域的R、G和B通道的平均值。
% 
% sum1 = zeros(1,3); 
% sum2 = zeros(1,3); 
% flag1 = zeros(1,3); 
% flag2 =zeros(1,3);
% for i = 1:m %遍历每一个像素点，可以根据需要选择自己需要处理的区域
%     for j = 1: n
%         if (hsv(i,j,3) < AverV) && (img(i,j,3) < AverB)
%             
%         elseif
%             
%         end
%         
%     end
% end



% 4.3. 分别计算每个通道的外部和内部的平均值的比率。
% 4.4. 将阴影区的每个像素与4.3中计算的常数相乘。
% 4.5. 对于阴影区的边缘像素，使用中值滤波来消除高频。

%% 
clear all, close all, clc

image = double(imread('images/56.jpg'))./255; % 颜色值归一化
s_im = size(image);

%*************************************************************************%

% 阴影检测
    % MASK: creating a shadow segmentation if no mask is available
    gray = rgb2gray(image);
    light_mask = double(bwareaopen(im2bw(gray, graythresh(gray)),200));
    h = fspecial('gaussian',20,0.5);
    light_mask = imfilter(light_mask,h);
    shadow_mask = 1 - light_mask;
    

    % SHADOW / LIGHT CORE DETECTION
    % structuring element
    struct_elem = [0 1 1 1 0; 1 1 1 1 1; 1 1 1 1 1; 1 1 1 1 1; 0 1 1 1 0];
    
    % shadow/light  core (morphology erode: pixels not on the blurred edge of the shadow area)
    shadow_core = imerode(shadow_mask, struct_elem);
    light_core = imerode(light_mask, struct_elem);
    % smoothing the mask
    smoothmask = conv2(shadow_mask, struct_elem/sum(sum(struct_elem)), 'same');
    
    % AVERAGE PIXEL INTENSITIES
    % shadow area
    shadow_avg_red = sum(sum(image(:,:,1).*shadow_core)) / sum(sum(shadow_core));
    shadow_avg_green = sum(sum(image(:,:,2).*shadow_core)) / sum(sum(shadow_core));
    shadow_avg_blue = sum(sum(image(:,:,3).*shadow_core)) / sum(sum(shadow_core));
    % light area
    light_avg_red = sum(sum(image(:,:,1).*light_core)) / sum(sum(light_core));
    light_avg_green = sum(sum(image(:,:,2).*light_core)) / sum(sum(light_core));
    light_avg_blue = sum(sum(image(:,:,3).*light_core)) / sum(sum(light_core));

% %     % K-MEANS CLUSTERING
%     im_hsv = rgb2hsv(image);
%     data = zeros(s_im(1)*s_im(2), 5);
%     for r=1:s_im(1)
%         for c=1:s_im(2)
%             data(r,:) = [im_hsv(r,c,1), im_hsv(r,c,2), im_hsv(r,c,3), r, c];
%         end
%     end
%     k = 5;
%     [IDX C] = kmeans(data, k, 'emptyaction', 'singleton');
%     [val, shadow_cluster] = min(C(:,2));
%     k_mask = vec2mat(IDX == shadow_cluster, s_im(2));
%     result_k_mask = im_hsv;
%     result_mask(:,:, 1) = im_hsv(:,:,1) .* k_mask;
%     result_mask(:,:, 2) = im_hsv(:,:,2) .* k_mask;
%     result_mask(:,:, 3) = im_hsv(:,:,3) .* k_mask;
%     figure, imshow(hsv2rgb(double(k_mask)))
