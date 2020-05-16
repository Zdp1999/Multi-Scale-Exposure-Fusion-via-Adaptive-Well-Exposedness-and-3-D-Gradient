clc, clear, close all;

imgSeqColor = loadImg('E:\codes\Kede_Ma\1\MEFDatabase\source image sequences\Venice_HDRsoft\'); % use im2double

tic
fI = spdmefStatic(imgSeqColor);
toc
imwrite(fI,'result_ma.png');
figure, imshow(fI);



