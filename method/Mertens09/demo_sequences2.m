clear all;close all;clc;
warning off;
addpath(genpath(pwd));

%%
folders_path = 'D:\MEF\dataset6\exposure_dataset_our\';
files = dir(folders_path);
length_file = length(files);
dataset6_Q = zeros(length_file - 2,1);
%%
for k = 3 : length_file
    images_path = [folders_path, files(k).name];
    imgs = load_images(images_path,1);
    I = imgs/255.0;
    img_result = exposure_fusion(I,[1 1 1]);
    want_to_save_result = 1;
    if (want_to_save_result)
        imwrite(img_result, ['E:\paper\dataset6_result\Mertens09\dataset6_result\',files(k).name,'_Mertens_CGF09.tif']);
    end
%% evalutate using MEF-SSIM
    imgSeqColor = uint8(load_images(images_path,1));
    [s1, s2, s3, s4] = size(imgSeqColor);
    imgSeq = zeros(s1, s2, s4);
    for n = 1:s4
        imgSeq(:, :, n) =  rgb2gray( squeeze( imgSeqColor(:,:,:,n) ) ); % color to gray conversion
    end
    fI1 = imread(['E:\paper\dataset6_result\Mertens09\dataset6_result\',files(k).name,'_Mertens_CGF09.tif']);
    fI1 = double(rgb2gray(fI1));
    [dataset6_Q(k-2,1), Qs1, QMap1] = mef_ms_ssim(imgSeq, fI1);   
end
save dataset6_Q.mat
