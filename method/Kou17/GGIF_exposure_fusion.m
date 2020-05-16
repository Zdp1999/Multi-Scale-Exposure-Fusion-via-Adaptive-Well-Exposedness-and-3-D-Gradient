function R = GGIF_exposure_fusion(I,m)
%   Code for
%   [1] "Multi-scale Exposure Fusion via Gradient Domain Guided Image Filtering", Fei Kou, Zhengguo Li, Changyun Wen, and Weihai Chen,  
%   in IEEE International Conference on Multimedia & Expo (ICME) 2017.

%   Code Author: Fei KOU
%   Email: koufei@hotmail.com
%   Date: 30/06/2017
%   Modified from the code provied by Dr. Tom Mertens et. al. for their paper 
%   [2] "Exposure Fusion" Tom Mertens, Jan Kautz and Frank Van Reeth, In Proceedings of Pacific Graphics 2007
%  
%   The code and the algorithm are for non-comercial use only.

% Usage:
%   result = GGIF_exposure_fusion(I,m);
%   Arguments:
%     'I': represents a stack of N color images (at double
%       precision). Dimensions are (height x width x 3 x N).
%     'm': 3-tuple that controls the per-pixel measures. The elements 
%     control contrast, saturation and well-exposedness, respectively.
%
BETA = 2; 
    
r = size(I,1);
c = size(I,2);
N = size(I,4);

nlev = floor(log(min(r,c)) / log(2))-BETA;
radius = 2^BETA; 
eps = 1/1024; 

%compute the measures and combines them into a weight map
contrast_parm = m(1);
sat_parm = m(2);
wexp_parm = m(3);

W = ones(r,c,N);
if (contrast_parm > 0)
    W = W.*contrast(I).^contrast_parm;
end
if (sat_parm > 0)
    W = W.*saturation(I).^sat_parm;
end
if (wexp_parm > 0)
    W = W.*well_exposedness(I).^wexp_parm;
end

%normalize weights: make sure that weights sum to one for each pixel
W = W + 1e-12;
W = W./repmat(sum(W,3),[1 1 N]);

Y = 0.299*I(:,:,1,:)+0.587*I(:,:,2,:)+0.114*I(:,:,3,:);
pyr = gaussian_pyramid(zeros(r,c,3),nlev);

for i = 1:N
    tmp =  Y(:,:,i);
    W(:,:,i) = W(:,:,i)*(mean(tmp(:))^2);
end
W = W./repmat(sum(W,3),[1 1 N]);

wsum = cell(nlev,1);
for i = 1:N
    pyrW = gaussian_pyramid(W(:,:,i),nlev);
    pyrY = gaussian_pyramid(Y(:,:,i),nlev);
    pyrI = laplacian_pyramid(I(:,:,:,i),nlev);
    
    for l = 1:nlev
        pyrW{l}=gguidedfilter(pyrY{l}, pyrW{l}, radius, eps);
    end
    
    for l = 1:nlev
        if i == 1
            wsum{l}= pyrW{l}+ 1e-12;
        else
            wsum{l} = wsum{l}+pyrW{l}+1e-12;
        end
        w = repmat(pyrW{l}+ 1e-12,[1 1 3]);
        pyr{l} = pyr{l} + w.*pyrI{l};
    end
end
for l = 1:nlev
    pyr{l} = pyr{l}./repmat(wsum{l},[1,1,3]);
end
R = reconstruct_laplacian_pyramid(pyr);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% contrast measure
function C = contrast(I)
h = [0 1 0; 1 -4 1; 0 1 0]; % laplacian filter
N = size(I,4);
C = zeros(size(I,1),size(I,2),N);
for i = 1:N
    mono = rgb2gray(I(:,:,:,i));
    C(:,:,i) = abs(imfilter(mono,h,'replicate'));
end

% saturation measure
function C = saturation(I)
N = size(I,4);
C = zeros(size(I,1),size(I,2),N);
for i = 1:N
    % saturation is computed as the standard deviation of the color channels
    R = I(:,:,1,i);
    G = I(:,:,2,i);
    B = I(:,:,3,i);
    mu = (R + G + B)/3;
    C(:,:,i) = sqrt(((R - mu).^2 + (G - mu).^2 + (B - mu).^2)/3);%./255;
end

% well-exposedness measure
function C = well_exposedness(I)
sig = .2; %%%only the color of final image will be effected by the value of sig.
N = size(I,4);
C = zeros(size(I,1),size(I,2),N);
for i = 1:N
    R = exp(-.5*(I(:,:,1,i) - .5).^2/sig.^2);
    G = exp(-.5*(I(:,:,2,i) - .5).^2/sig.^2);
    B = exp(-.5*(I(:,:,3,i) - .5).^2/sig.^2);
    C(:,:,i) = R.*G.*B;
end
