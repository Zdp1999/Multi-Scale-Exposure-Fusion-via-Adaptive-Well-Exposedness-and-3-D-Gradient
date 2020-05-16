clear all;close all;clc;
warning off;
file_dir = fileparts(mfilename('fullpath')); cd(file_dir);
addpath(genpath(pwd));
dataset1_1_Q = zeros(360,1);

%% add path of a folder including support functions
%addpath('functions');

%% read multi-exposed rgb image sequence (scaling to [0,1])
for i = 1 : 360
    i
    image_name = num2str(i);
    imgs_rgb = load_images(['D:\MEF\MEF_dataset1\', image_name, '\'],1);
    I = imgs_rgb/255.0;
%tic;
%% compute luminance image of rgb image sequence
img_result = GGIF_exposure_fusion(I,[1 1 1]); 

%% show and save result image (optional)
% figure();
% imshow(img_result);
%toc;
want_to_save_result = 1;
if (want_to_save_result)
    imwrite(img_result, ['E:\paper\dataset1_result\Kou17\Dataset_Part1\',image_name,'_kou_ICME17.tif']);
    %imwrite(img_result, ['G:/dataset/HDR/MEFDatabase/result/',image_name,'/result01.png']);
end

%% evalutate using MEF-SSIM
    imgSeqColor = uint8(load_images(['D:\MEF\MEF_dataset1\', image_name, '\'],1));
    [s1, s2, s3, s4] = size(imgSeqColor);
    imgSeq = zeros(s1, s2, s4);
    for n = 1:s4
        imgSeq(:, :, n) =  rgb2gray( squeeze( imgSeqColor(:,:,:,n) ) ); % color to gray conversion
    end
    fI1 = imread(['E:\paper\dataset1_result\Kou17\Dataset_Part1\',image_name,'_kou_ICME17.tif']);
    fI1 = double(rgb2gray(fI1));
    [dataset1_1_Q(i,1), Qs1, QMap1] = mef_ms_ssim(imgSeq, fI1);   
end
save Q1.mat

%%%
clear all;close all;clc;
warning off;
file_dir = fileparts(mfilename('fullpath')); cd(file_dir);
addpath(genpath(pwd));
Q2 = zeros(229,1);

%% add path of a folder including support functions
%addpath('functions');

%% read multi-exposed rgb image sequence (scaling to [0,1])
for i = 1 : 229
    image_name = num2str(i);
    imgs_rgb = load_images(['D:\MEF\MEF_dataset2\', image_name, '\'],1);
    I = imgs_rgb/255.0;
% tic;
%% compute luminance image of rgb image sequence
img_result = GGIF_exposure_fusion(I,[1 1 1]); 

%% show and save result image (optional)
% figure();
% imshow(img_result);
% toc;
want_to_save_result = 1;
if (want_to_save_result)
    imwrite(img_result, ['E:\paper\dataset1_result\Kou17\Dataset_Part2\',image_name,'_kou_ICME17.tif']);
    %imwrite(img_result, ['G:/dataset/HDR/MEFDatabase/result/',image_name,'/result01.png']);
end

%% evalutate using MEF-SSIM
    imgSeqColor = uint8(load_images(['D:\MEF\MEF_dataset2\', image_name, '\'],1));
    [s1, s2, s3, s4] = size(imgSeqColor);
    imgSeq = zeros(s1, s2, s4);
    for n = 1:s4
        imgSeq(:, :, n) =  rgb2gray( squeeze( imgSeqColor(:,:,:,n) ) ); % color to gray conversion
    end
    fI1 = imread(['E:\paper\dataset1_result\Kou17\Dataset_Part2\',image_name,'_kou_ICME17.tif']);
    fI1 = double(rgb2gray(fI1));
    [Q2(i,1), Qs1, QMap1] = mef_ms_ssim(imgSeq, fI1);   
end
save Q2.mat