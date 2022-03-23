% im_ycbcr COLORSPACE:
% Hybridation of the additive and light model based methods on im_ycbcr colourspace.
%

clear all, close all, clc

image = double(imread('images/redtest.png'))./255; % 颜色值归一化
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

% Method 4: COMBINED ADDITIVE AND LIGHT MODEL BASED SHADOW REMOVAL IN im_ycbcr COLOURSPACE
% conversion to YCbCr colorspace
im_ycbcr = rgb2ycbcr(image);
% computing averade channel values in im_ycbcr space
shadow_avg_y = sum(sum(im_ycbcr(:,:,1).*shadow_core)) / sum(sum(shadow_core));
shadow_avg_cb = sum(sum(im_ycbcr(:,:,2).*shadow_core)) / sum(sum(shadow_core));
shadow_avg_cr = sum(sum(im_ycbcr(:,:,3).*shadow_core)) / sum(sum(shadow_core));
%
litavg_y = sum(sum(im_ycbcr(:,:,1).*light_core)) / sum(sum(light_core));
litavg_cb = sum(sum(im_ycbcr(:,:,2).*light_core)) / sum(sum(light_core));
litavg_cr = sum(sum(im_ycbcr(:,:,3).*light_core)) / sum(sum(light_core));
% computing ratio, and difference in im_ycbcr space
diff_y = litavg_y - shadow_avg_y;
diff_cb = litavg_cb - shadow_avg_cb;
diff_cr = litavg_cr - shadow_avg_cr;

ratio_y = litavg_y/shadow_avg_y;
ratio_cb = litavg_cb/shadow_avg_cb;
ratio_cr = litavg_cr/shadow_avg_cr;
% shadow correction: Y->additive, Cb&Cr-> basic light model
aux_result_im_ycbcr = im_ycbcr;
aux_result_im_ycbcr(:,:,1) = im_ycbcr(:,:,1) + shadow_mask * diff_y;
aux_result_im_ycbcr(:,:,2) = im_ycbcr(:,:,2).*light_mask + shadow_mask.*ratio_cb.*im_ycbcr(:,:,2);
aux_result_im_ycbcr(:,:,3) = im_ycbcr(:,:,3).*light_mask + shadow_mask.*ratio_cr.*im_ycbcr(:,:,3);
% conversion back to rgb colourspace
result_im_ycbcr = ycbcr2rgb(aux_result_im_ycbcr);

imshow(smoothmask), title('Smooth Mask')

imshow(result_im_ycbcr), title('method 2')