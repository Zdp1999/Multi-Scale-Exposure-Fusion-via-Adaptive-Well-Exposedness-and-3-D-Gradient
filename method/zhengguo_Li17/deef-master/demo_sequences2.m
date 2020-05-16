clear all;close all;clc;
warning off;
addpath(genpath(pwd));

%%
folders_path = 'D:\MEF\dataset6\exposure_dataset_our\';
files = dir(folders_path);
length_file = length(files);
dataset6_Q = zeros(length_file - 2,1);
%%
for N = 3 : length_file
     images_path = [folders_path, files(N).name];
%     imgs = load_images(images_path,1);
%     nlev = 2;
%     img_result = exposure_fusion(imgs,nlev);
    [filenames, num_LDR] = readDir2([images_path, '\']);
    
    temp = imread(filenames{1});
    [height,width,color]=size(temp);
    
    I = zeros(height,width,color,num_LDR);
    for k = 1:num_LDR
        temp = double(imread(filenames{k}));
        for i = 1:height
            for j = 1:width
                for n = 1:color
                    I(i,j,n,k)=temp(i,j,n)/255;
                end
            end
        end
    end
    
    clear temp;
    
    %tic
    %%%WGIF based exposure fusion
    R = exposure_fusion(I,[1 1 1]); %% [1,1,1] to enable the three weighting factors.
    % figure;
    % imshow(R);
    % imwrite(R, 'fusion_base.png','png');
    % %%%Detail enhancement
    for i = 1:height
        for j = 1:width
            for n = 1:color
                for k = 1:num_LDR
                    I(i,j,n,k) = floor(I(i,j,n,k)*255);
                end
                R(i,j,n) = max(min(255, floor(R(i,j,n)*255)),0);
            end
        end
    end
    R_d = extract_details(I, R);
    theta = 1;
    Y_d = 2.^(theta*R_d);
    
    for k=1:1:3
        R(:,:,k) = floor(R(:,:,k).*(Y_d));
    end
    img_result = R/255;
    want_to_save_result = 1;
    if (want_to_save_result)
        imwrite(img_result, ['E:\paper\dataset6_result\zhengguo_Li17\dataset6_result\',files(N).name,'_zhengguo_Li_TIP17.tif']);
    end
%% evalutate using MEF-SSIM
    imgSeqColor = uint8(load_images1(images_path,1));
    [s1, s2, s3, s4] = size(imgSeqColor);
    imgSeq = zeros(s1, s2, s4);
    for n = 1:s4
        imgSeq(:, :, n) =  rgb2gray( squeeze( imgSeqColor(:,:,:,n) ) ); % color to gray conversion
    end
    fI1 = imread(['E:\paper\dataset6_result\zhengguo_Li17\dataset6_result\',files(N).name,'_zhengguo_Li_TIP17.tif']);
    fI1 = double(rgb2gray(fI1));
    [dataset6_Q(N-2,1), Qs1, QMap1] = mef_ms_ssim(imgSeq, fI1);   
end
save dataset6_Q.mat
