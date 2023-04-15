% BASIC LIGHT MODEL:
% Basic light model containing an ambient & a directed light ratio.

clear all, close all, clc

% [filename,pathname]=uigetfile({'*.jpg;*.bmp;*.tif;*.png;*.gif','All Image Files';'*.*','All Files'});
% image = imread([pathname,filename]);


% image = double(imread('images/redbig.png'))./255; % 颜色值归一化
image = double(imread('images/56.jpg'))./255;
s_im = size(image);

%*************************************************************************%

% 阴影检测
    % MASK: creating a shadow segmentation if no mask is available
    gray = rgb2gray(image);
    light_mask = double(bwareaopen(im2bw(gray, graythresh(gray)),200));
%     h = fspecial('gaussian',20,0.5);    
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

%     % K-MEANS CLUSTERING
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


% Method 2: BASIC , LIGHT MODEL BASED SHADOW REMOVAL
result_basic_model = zeros(s_im);
% computing ratio of shadow/lit area luminance
ratio_red = light_avg_red/shadow_avg_red;
ratio_green = light_avg_green/shadow_avg_green;
ratio_blue = light_avg_blue/shadow_avg_blue;
%
result_basic_model(:,:,1) = (light_mask + shadow_mask.*ratio_red).*image(:,:,1);
result_basic_model(:,:,2) = image(:,:,2).*light_mask + shadow_mask.*ratio_green.*image(:,:,2);
result_basic_model(:,:,3) = image(:,:,3).*light_mask + shadow_mask.*ratio_blue.*image(:,:,3);

imshow(shadow_mask), title('Shadow Mask')
figure

imshow(result_basic_model), title('Basic')

% imwrite(result_basic_model,'methodBasic.png')
imwrite(result_basic_model,'basic.png')