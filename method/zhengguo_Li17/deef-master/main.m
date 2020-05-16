%function example;
clear all;
close all;

tic;
dirName = ('tower\');
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

% figure;
% imshow(R/255);
% 
% runtime = toc;
% imwrite(R/255, 'fusion_enhanced.png','png');
% fprintf('the running time is %f\n', runtime);


