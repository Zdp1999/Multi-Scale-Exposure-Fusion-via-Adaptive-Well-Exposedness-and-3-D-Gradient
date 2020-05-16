clear all;close all;clc;
warning off;
file_dir = fileparts(mfilename('fullpath')); cd(file_dir);
addpath(genpath(pwd));
Q1 = zeros(360,1);

%% add path of a folder including support functions
%addpath('functions');

%% read multi-exposed rgb image sequence (scaling to [0,1])
for N = 1 : 360
    N
    image_name = num2str(N);
    %imgs_rgb = load_images(['D:\MEF\MEF_dataset1\', image_name, '\'],1);
    dirName = (['D:\data\MEF_dataset1\', image_name, '\']);
    [filenames, num_LDR] = readDir2(dirName);
    
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
        imwrite(img_result, ['D:\paper\zhengguo_Li17\Dataset_Part1\',image_name,'_zhengguo_Li_TIP17.tif']);
        %imwrite(img_result, ['G:/dataset/HDR/MEFDatabase/result/',image_name,'/result01.png']);
    end

%% evalutate using MEF-SSIM
    imgSeqColor = uint8(load_images(['D:\data\MEF_dataset1\', image_name, '\'],1));
    [s1, s2, s3, s4] = size(imgSeqColor);
    imgSeq = zeros(s1, s2, s4);
    for n = 1:s4
        imgSeq(:, :, n) =  rgb2gray( squeeze( imgSeqColor(:,:,:,n) ) ); % color to gray conversion
    end
    fI1 = imread(['D:\paper\zhengguo_Li17\Dataset_Part1\',image_name,'_zhengguo_Li_TIP17.tif']);
    fI1 = double(rgb2gray(fI1));
    [Q1(N,1), Qs1, QMap1] = mef_ms_ssim(imgSeq, fI1);   
end
save Q1.mat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;close all;clc;
warning off;
file_dir = fileparts(mfilename('fullpath')); cd(file_dir);
addpath(genpath(pwd));
Q2 = zeros(229,1);

%% add path of a folder including support functions
%addpath('functions');

%% read multi-exposed rgb image sequence (scaling to [0,1])
for N = 1 : 229
    N
    image_name = num2str(N);
    %imgs_rgb = load_images(['D:\MEF\MEF_dataset1\', image_name, '\'],1);
%     dirName = (['D:\data\MEF_dataset2\', image_name, '\']);
%     [filenames, num_LDR] = readDir2(dirName);
%     
%     temp = imread(filenames{1});
%     [height,width,color]=size(temp);
%     
%     I = zeros(height,width,color,num_LDR);
%     for k = 1:num_LDR
%         temp = double(imread(filenames{k}));
%         for i = 1:height
%             for j = 1:width
%                 for n = 1:color
%                     I(i,j,n,k)=temp(i,j,n)/255;
%                 end
%             end
%         end
%     end
%     
%     clear temp;
%     
%     %tic
%     %%%WGIF based exposure fusion
%     R = exposure_fusion(I,[1 1 1]); %% [1,1,1] to enable the three weighting factors.
%     % figure;
%     % imshow(R);
%     % imwrite(R, 'fusion_base.png','png');
%     % %%%Detail enhancement
%     for i = 1:height
%         for j = 1:width
%             for n = 1:color
%                 for k = 1:num_LDR
%                     I(i,j,n,k) = floor(I(i,j,n,k)*255);
%                 end
%                 R(i,j,n) = max(min(255, floor(R(i,j,n)*255)),0);
%             end
%         end
%     end
%     R_d = extract_details(I, R);
%     theta = 1;
%     Y_d = 2.^(theta*R_d);
%     
%     for k=1:1:3
%         R(:,:,k) = floor(R(:,:,k).*(Y_d));
%     end
%     img_result = R/255;
% 
%     want_to_save_result = 1;
%     if (want_to_save_result)
%         imwrite(img_result, ['D:\paper\zhengguo_Li17\Dataset_Part2\',image_name,'_zhengguo_Li_TIP17.tif']);
%         %imwrite(img_result, ['G:/dataset/HDR/MEFDatabase/result/',image_name,'/result01.png']);
%     end

%% evalutate using MEF-SSIM
    imgSeqColor = uint8(load_images(['D:\data\MEF_dataset2\', image_name, '\'],1));
    [s1, s2, s3, s4] = size(imgSeqColor);
    imgSeq = zeros(s1, s2, s4);
    for n = 1:s4
        imgSeq(:, :, n) =  rgb2gray( squeeze( imgSeqColor(:,:,:,n) ) ); % color to gray conversion
    end
    fI1 = imread(['D:\paper\zhengguo_Li17\Dataset_Part2\',image_name,'_zhengguo_Li_TIP17.tif']);
    fI1 = double(rgb2gray(fI1));
    [Q2(N,1), Qs1, QMap1] = mef_ms_ssim(imgSeq, fI1);   
end
save Q2.mat

