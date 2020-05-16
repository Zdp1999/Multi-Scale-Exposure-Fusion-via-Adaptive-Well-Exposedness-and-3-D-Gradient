clear all;close all;

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

%% read multi-exposed rgb image sequence (scaling to [0,1])
T = 0;
for i = 12%1 : 21
    image_name = image_names{i};
    imgs_rgb = load_images(['E:\second_paper\source_images\', image_name],1);
    %imgs_rgb = imgs_rgb/255.0;
    tic;
    img_result = GFF(imgs_rgb);

%% show and save result image (optional)
% figure();
% imshow(img_result);
toc;
want_to_save_result = 1;
if (want_to_save_result)
    imwrite(img_result, ['E:\second_paper\result\4\',image_name,'_shutao_13.tif']);
    %imwrite(img_result, ['G:/dataset/HDR/MEFDatabase/result/',image_name,'/result01.png']);
end

%% evalutate using MEF-SSIM
    imgSeqColor = uint8(load_images(['E:\second_paper\source_images\', image_name],1));
    [s1, s2, s3, s4] = size(imgSeqColor);
    imgSeq = zeros(s1, s2, s4);
    for n = 1:s4
        imgSeq(:, :, n) =  rgb2gray( squeeze( imgSeqColor(:,:,:,n) ) ); % color to gray conversion
    end
    fI1 = imread(['E:\second_paper\result\4\',image_name,'_shutao_13.tif']);
    fI1 = double(rgb2gray(fI1));
    [Q(i,1), Qs1, QMap1] = mef_ms_ssim(imgSeq, fI1);   
end
xlswrite('E:\second_paper\result\4\score.xlsx',Q)