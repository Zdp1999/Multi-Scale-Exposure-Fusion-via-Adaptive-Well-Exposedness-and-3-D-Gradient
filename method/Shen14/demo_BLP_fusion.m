% function demo_BLP_fusion

clear all;close all;clc;
addpath('code');
I= load_images('images\Set');
% input_folder = 'E:\codes\Kede_Ma\3\database_release\source_images\Set\';
% I= load_images(input_folder);
%%
% recommend small size for testing
% I = min(1,max(0, imresize(I,1/16) )); 
% fprintf('processing %d images sequence:land\n',size(I,4));
% nlev: pyramid level
nlev=2;
tic
R = exposure_fusion(I,nlev);
toc
% figure('Name','Result'); 
% imshow(R); 
imwrite(R,'result.png');

