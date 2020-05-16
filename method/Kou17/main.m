clear all;
close all;

dirName = ('treeunil\');
delete ('treeunil\thumbs.db');
[filenames, num_LDR] = readDir2(dirName);
temp = imread(filenames{1});
[height,width,color]=size(temp);
I = zeros(height,width,color,num_LDR);
for k = 1:num_LDR
    temp = double(imread(filenames{k}));
    I(:,:,:,k)=temp/255;
end

clear temp;

%exposure fusion by the method in ¡°Multi-scale Exposure Fusion via Gradient Domain Guided Image Filtering¡±
R_GGIF = GGIF_exposure_fusion(I,[1 1 1]); 
figure('Name','FusionResult_ICME2017');
imshow(R_GGIF);
imwrite(R_GGIF, 'Treeunil_ICME17.png','png');

% %exposure fusion by T. Mertens's method
% R = exposure_fusion(I,[1 1 1]); 
% figure('Name','FusionResult_T.Mertens');
% imshow(R);
% imwrite(R, 'Treeunil_Mertens.png','png');


