% ENHANCED LIGHT MODEL:
% Enhanced shadow removal model containing an ambient & a directed light ratio
% and shadow coefficient.

clear all, close all, clc

image = double(imread('images/redbig.png'))./255; % 颜色值归一化
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



% Method : ADVANCE, LIGHT MODEL BASED SHADOW REMOVAL
result_enhanced_model = zeros(s_im);
% computing ratio of the luminances of the directed, and global lights
ratio_red = light_avg_red/shadow_avg_red - 1;
ratio_green = light_avg_green/shadow_avg_green - 1;
ratio_blue = light_avg_blue/shadow_avg_blue - 1;
% applying shadow removal
result_enhanced_model(:,:,1) = (ratio_red + 1)./((1-smoothmask)*ratio_red + 1).*image(:,:,1);
result_enhanced_model(:,:,2) = (ratio_green + 1)./((1-smoothmask)*ratio_green + 1).*image(:,:,2);
result_enhanced_model(:,:,3) = (ratio_blue + 1)./((1-smoothmask)*ratio_blue + 1).*image(:,:,3);

imshow(shadow_mask), title('Shadow Mask');
figure

imshow(result_enhanced_model), title('Light')

imwrite(result_enhanced_model,'methodLight.png')