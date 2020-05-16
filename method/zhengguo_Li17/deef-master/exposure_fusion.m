%
% Implementation of Exposure Fusion
%
% written by Tom Mertens, Hasselt University, August 2007
% e-mail: tom.mertens@gmail.com
%
% This work is described in
%   "Exposure Fusion"
%   Tom Mertens, Jan Kautz and Frank Van Reeth
%   In Proceedings of Pacific Graphics 2007
%
%
% Usage:
%   result = exposure_fusion(I,m);
%   Arguments:
%     'I': represents a stack of N color images (at double
%       precision). Dimensions are (height x width x 3 x N).
%     'm': 3-tuple that controls the per-pixel measures. The elements 
%     control contrast, saturation and well-exposedness, respectively.
%
% Example:
%   'figure; imshow(exposure_fusion(I, [0 0 1]);'
%   This displays the fusion of the images in 'I' using only the well-exposedness
%   measure
%

function R = exposure_fusion(I,m)

method = 2; %%%1: the Gaussian pyramid; 2: the WGIF based pyramid; 3: the hybrid pyramid
Kappa = 3;
BETA = 2;%2; %%%test four different choices of BETA as 0, 1, 2, 3
% Alpha_l = 255; %127;%%control weight of the darkest image

    
r = size(I,1);
c = size(I,2);
N = size(I,4);

nlev = floor(log(min(r,c)) / log(2))-BETA;
Kappa = min(Kappa, nlev-1);
radius = 2^BETA; 
eps = 1/1024;    

W = ones(r,c,N);
Y = zeros(r,c,N);

 %compute the measures and combines them into a weight map
contrast_parm = m(1);
sat_parm = m(2);
wexp_parm = m(3);

% W_Lo = ones(256,1);
% 
% 
% 
% %%%For the darkest image. Use the Y component
% for z=1:256
%     W_Lo(z) = 1+Alpha_l*(z/256)^0.25;
% end


if (contrast_parm > 0)
    W = W.*contrast(I).^contrast_parm;
end
if (sat_parm > 0)
    W = W.*saturation(I).^sat_parm;
end
if (wexp_parm > 0)
    W = W.*well_exposedness(I).^wexp_parm;
end
% Y_ave = zeros(N,1);    
for i=1:N
    Y(:,:,i) = 0.299*I(:,:,1,i)+0.587*I(:,:,2,i)+0.114*I(:,:,3,i);
%     Y_ave(i) = mean(mean(Y(:,:,i)))^POW_T;
%     W(:,:,i) = W(:,:,i)*Y_ave(i);
%     if i==1
%         for ii=1:r
%             for jj=1:c
%                 W(ii,jj,i) = W(ii,jj,i)*W_Lo(floor(Y(ii,jj,i)+0.5)+1);
%             end
%         end
%     end        
end
%normalize weights: make sure that weights sum to one for each pixel
W = W + 1e-12;
W = W./repmat(sum(W,3),[1 1 N]);
if method==1 %%%Gaussian pyramid
    % create empty pyramid
    nlev = floor(log(min(r,c)) / log(2));
    pyr = gaussian_pyramid(zeros(r,c,3),nlev);
    
    % multiresolution blending
    for i = 1:N
        % construct pyramid from each input image
        pyrW = gaussian_pyramid(W(:,:,i),nlev);
        pyrI = laplacian_pyramid(I(:,:,:,i),nlev);
        
        % blend
        for l = 1:nlev
            w = repmat(pyrW{l},[1 1 3]);
            pyr{l} = pyr{l} + w.*pyrI{l};
        end
    end
    
    % reconstruct
    R = reconstruct_laplacian_pyramid(pyr);
elseif method==2 %%%Fuse all color components via simplified WGIF pyramids
    % create empty pyramid
    pyr = gaussian_pyramid(zeros(r,c,3),nlev);
    wsum = cell(nlev,1);
    odd = zeros(2,1);
    filter = pyramid_filter;
    % multiresolution blending
    for i = 1:N
        % construct pyramid from each input image
        pyrW = gaussian_pyramid(W(:,:,i),nlev);
        [pyrY, pyr_Ref] = laplacian_pyramid_layer(Y(:,:,i),nlev);
        pyrI = laplacian_pyramid(I(:,:,:,i),nlev);
        [pyrW{nlev}, a, b] = guidedfilter_WMSE_FixedRadius(pyr_Ref{nlev}, pyrW{nlev}, radius, eps);        
        for l = nlev-1: -1: Kappa+2
            [hei_t_s, wid_t_s] = size(pyrW{l+1}); 
            [hei_t, wid_t] = size(pyrW{l});
            odd(1) = 2*hei_t_s-hei_t;
            odd(2) = 2*wid_t_s-wid_t;
            a = upsample(a,odd,filter);
            b = upsample(b,odd,filter);           
            pyrW{l} = a.*pyr_Ref{l}+b;
            w_tmp = pyrW{l}+1e-12;
            if i == 1 
                wsum{l} = w_tmp;
            else
                wsum{l} = wsum{l}+w_tmp;
            end            
            w = repmat(w_tmp,[1 1 3]);            
            pyr{l} = pyr{l} + w.*pyrI{l};
        end
        clear a;
        clear b;
        [pyrW{Kappa+1}, a, b] = guidedfilter_WMSE_FixedRadius(pyr_Ref{Kappa+1}, pyrW{Kappa+1}, radius, eps); 
        for l = Kappa: -1: 1
            [hei_t_s, wid_t_s] = size(pyrW{l+1}); 
            [hei_t, wid_t] = size(pyrW{l});
            odd(1) = 2*hei_t_s-hei_t;
            odd(2) = 2*wid_t_s-wid_t;
            a = upsample(a,odd,filter);
            b = upsample(b,odd,filter);           
            pyrW{l} = a.*pyr_Ref{l}+b;       
            w_tmp = pyrW{l}+1e-12;
            if i == 1 
                wsum{l}= w_tmp;
            else
                wsum{l} = wsum{l}+w_tmp;
            end   
            w = repmat(w_tmp,[1 1 3]); 
            pyr{l} = pyr{l} + w.*pyrI{l};
        end        
        if i==1
            w_tmp = pyrW{nlev}+ 1e-12;
            wsum{nlev} = w_tmp;
            w = repmat(w_tmp,[1 1 3]); 
            pyr{nlev} = pyr{nlev} + w.*pyrI{nlev};
            w_tmp = pyrW{Kappa+1}+ 1e-12;
            wsum{Kappa+1} = w_tmp;
            w = repmat(w_tmp,[1 1 3]);
            pyr{Kappa+1} = pyr{Kappa+1} + w.*pyrI{Kappa+1};            
        else
            w_tmp = pyrW{nlev}+ 1e-12;
            wsum{nlev} = wsum{nlev}+w_tmp;
            w = repmat(w_tmp,[1 1 3]);
            pyr{nlev} = pyr{nlev} + w.*pyrI{nlev}; 
            w_tmp = pyrW{Kappa+1}+ 1e-12;
            wsum{Kappa+1} = wsum{Kappa+1}+w_tmp;
            w = repmat(w_tmp,[1 1 3]);
            pyr{Kappa+1} = pyr{Kappa+1} + w.*pyrI{Kappa+1};          
        end  
        
    end
    for l = 1: nlev
        w = repmat(wsum{l},[1 1 3]);
        pyr{l} = pyr{l}./w;
    end

    
    % reconstruct
    R = reconstruct_laplacian_pyramid(pyr);
else %%%hybrid pyramid
 %%%Fuse all color components via simplified WGIF pyramids
    % create empty pyramid
    pyr = gaussian_pyramid(zeros(r,c,3),nlev);
    wsum = cell(nlev,1);
    odd = zeros(2,1);
    filter = pyramid_filter;
    LAYER = 2;
    % multiresolution blending
    for i = 1:N
        % construct pyramid from each input image
        pyrW = gaussian_pyramid(W(:,:,i),nlev);
        [pyrY, pyr_Ref] = laplacian_pyramid_layer(Y(:,:,i),nlev);
        pyrI = laplacian_pyramid(I(:,:,:,i),nlev);
        [pyrW{nlev}, a, b] = guidedfilter_WMSE_FixedRadius(pyr_Ref{nlev}, pyrW{nlev}, radius, eps);
        for l = nlev-1: -1: nlev-LAYER
            [hei_t_s, wid_t_s] = size(pyrW{l+1}); 
            [hei_t, wid_t] = size(pyrW{l});
            odd(1) = 2*hei_t_s-hei_t;
            odd(2) = 2*wid_t_s-wid_t;
            a = upsample(a,odd,filter);
            b = upsample(b,odd,filter);           
            pyrW{l} = a.*pyr_Ref{l}+b;
        end        
        clear a;
        clear b; 
        for l = nlev: -1: 1
            w_tmp = pyrW{l}+1e-12;
            if i == 1 
                wsum{l} = w_tmp;
            else
                wsum{l} = wsum{l}+w_tmp;
            end            
            w = repmat(w_tmp,[1 1 3]);            
            pyr{l} = pyr{l} + w.*pyrI{l};
        end
    end
    for l = 1: nlev
        w = repmat(wsum{l},[1 1 3]);
        pyr{l} = pyr{l}./w;
    end
    % reconstruct
    R = reconstruct_laplacian_pyramid(pyr);
end




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
