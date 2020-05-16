clear all;close all;
clc;
image_names = cell(21,1);% 21¸ösequences
Q = zeros(21,1);
image_names{1} = 'Arno';
image_names{2} = 'Balloons';
image_names{3} = 'BelgiumHouse';
image_names{4} = 'Cave';
image_names{5} = 'ChineseGarden';
image_names{6} = 'Church';
image_names{7} = 'Farmhouse';
image_names{8} = 'House';
image_names{9} = 'Lamp';
image_names{10} = 'Landscape';
image_names{11} = 'Laurenziana';
image_names{12} = 'MadisonCapitol';
image_names{13} = 'Mask';
image_names{14} = 'Office';
image_names{15} = 'Ostrow';
image_names{16} = 'Room';
image_names{17} = 'Set';
image_names{18} = 'Tower';
image_names{19} = 'Venice';
image_names{20} = 'Window';
image_names{21} = 'YellowHall';


%% add path of a folder including support functions
addpath('functions');
%% read multi-exposed rgb image sequence (scaling to [0,1])
for i = 12%1 : 21
    image_name = image_names{i};
    I = load_images(['E:\second_paper\source_images\',image_name],1);
    I = I/255.0;
    %imgs_rgb = imgs_rgb/255.0;
    tic;
%% compute luminance image of rgb image sequence
% imgs_lum = rgb2lum(imgs_rgb);
% %% ×î¼Ñ
% [w2, lev] = weight_8(imgs_lum);
% %% compute weight2 using luminance gradient
% w3 = weight_2(imgs_rgb);
% %% corporate weight1 & weight2 and refine weight with wlsFilter
% p1 = 1;
% p2 = 1;
% p3 = 2; % 2
% %w = (w1*0.4 + w4*0.6) .* w3;
% %w = (w1.^p1).*(w3.^p3);
% w = (w2.^p2).*(w3.^p3);
% % w = feature_map(w);
% %%
% w = refine_weight(w);
% % w = refine(w, imgs_lum);
% %% fuse images using pyramid decomposition
% %lev = 7; % 7
% img_result = fusion_pyramid(imgs_rgb, w, lev);
R = exposure_fusion(I,[1 1 1]);
%% show and save result image (optional)
% figure();
% imshow(img_result);
toc;
want_to_save_result = 1;
if (want_to_save_result)
    imwrite(R, ['E:\second_paper\result\1\',image_name,'_Mertens.tif']);
    %imwrite(img_result, ['G:/dataset/HDR/MEFDatabase/result/',image_name,'/result01.png']);
end

%% evalutate using MEF-SSIM
    imgSeqColor = uint8(load_images(['E:\second_paper\source_images\', image_name],1));
    [s1, s2, s3, s4] = size(imgSeqColor);
    imgSeq = zeros(s1, s2, s4);
    for n = 1:s4
        imgSeq(:, :, n) =  rgb2gray( squeeze( imgSeqColor(:,:,:,n) ) ); % color to gray conversion
    end
    fI1 = imread(['E:\second_paper\result\1\',image_name,'_Mertens.tif']);
    fI1 = double(rgb2gray(fI1));
    [Q(i,1), Qs1, QMap1] = mef_ms_ssim(imgSeq, fI1);   
end
%xlswrite('E:\second_paper\result\1\score.xlsx',Q)