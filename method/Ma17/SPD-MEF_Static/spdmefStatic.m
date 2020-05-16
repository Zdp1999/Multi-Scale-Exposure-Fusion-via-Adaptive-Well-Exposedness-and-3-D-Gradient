function fI = spdmefStatic(imgSeqColor, varargin)
% ========================================================================
% Structural patch decomposition based multi-exposure image fusion (SPD-MEF)
% Static Version 1.0
% Copyright(c) 2015 Kede Ma and Zhou Wang
% All Rights Reserved.
%
% For an advanced version that handles camera and object motion, please
% refer to the following paper:
% K. Ma et al., "Robust Multi-Exposure Image Fusion:
% A Structural Patch Decomposition Approach," submitted to 
% IEEE Transactions on Image Processing.
% ----------------------------------------------------------------------
% Permission to use, copy, or modify this software and its documentation
% for educational and research purposes only and without fee is hereby
% granted, provided that this copyright notice and the original authors'
% names appear on all copies and supporting documentation. This program
% shall not be used, rewritten, or adapted as the basis of a commercial
% software or hardware product without first obtaining permission of the
% authors. The authors make no representations about the suitability of
% this software for any purpose. It is provided "as is" without express
% or implied warranty.
%----------------------------------------------------------------------
% This is an implementation of SPD-MEF for producing static MEF images.
% Please refer to the following paper:
%
% K. Ma and Z. Wang, "Multi-Exposure Image Fusion: A Patch-wise Approach," 
% in IEEE International Conference on Image Processing (ICIP), 2015.
%
% Kindly report any suggestions or corrections to k29ma@uwaterloo.ca
%
%----------------------------------------------------------------------
%
% Required Input : (1) imgSeqColor: source image sequence in [0-1] RGB.
% Optional Input : (2) p: the exponent parameter.  default value p = 4;
%                  (3) gSig: the spread of the global Gaussian. 
%                      defualt value: gSig = 0.2;
%                  (4) lSig: the spread of the local Gaussian. 
%                      defualt value: lSig = 0.5;
%                  (5) wSize: the patch size. defualt value: wSize = 11;
%                  (6) stepSize: the stride. 
%                      defualt value: stepSize = floor(wSize/10) = 2;
% Output:          (1) fI: The fused image.
%
% Basic Usage:
%   Given the input source image sequence imgColorSeq
%
%   fI = spdmefStatic(imgColorSeq);
%
% Advanced Usage:
%   
% fI = spdmefStatic(imgColorSeq, 'wSize', 9);
% 
%========================================================================
% input params parsing
    params = inputParser;
    
    default_p = 4;
    default_gSig = 0.2;
    default_lSig = 0.5;
    default_wSize = 11;
    default_stepSize = 2;
    
    addRequired(params,'imgSeqColor');
    addParameter(params, 'p', default_p, @isnumeric);
    addParameter(params, 'gSig', default_gSig, @isnumeric);
    addParameter(params, 'lSig', default_lSig, @isnumeric);
    addParameter(params, 'wSize', default_wSize, @isnumeric);
    addParameter(params, 'stepSize', default_stepSize, @isnumeric);
    
    parse(params, imgSeqColor, varargin{:});
    
% initialization 
    p = params.Results.p;
    gSig = params.Results.gSig;
    lSig = params.Results.lSig;
    wSize = params.Results.wSize;
    stepSize = params.Results.stepSize;
    
    window = ones(wSize);
    window = window / sum(window(:));

    imgSeqColor = double(imgSeqColor);
    [s1, s2, s3, s4] = size(imgSeqColor); 
    xIdxMax = s1-wSize+1;
    yIdxMax = s2-wSize+1;
    
    
% computing statistics
    gMu = zeros(xIdxMax, yIdxMax, s4); % global mean intensity
    for i = 1 : s4
        img = imgSeqColor(:,:,:,i);
        gMu(:,:,i) = ones(xIdxMax, yIdxMax) * mean(img(:));
    end
   
    temp  = zeros(xIdxMax, yIdxMax, s3);
    lMu   = zeros(xIdxMax, yIdxMax, s4); % local mean intensity
    lMuSq = zeros(xIdxMax, yIdxMax, s4);
    for i = 1 : s4
        for j = 1 : s3
            temp(:,:,j) = filter2(window, imgSeqColor(:, :, j, i), 'valid');
        end
        lMu(:,:,i) = mean(temp, 3); % (R + G + B) / 3;
        lMuSq(:,:,i) = lMu(:,:,i) .* lMu(:,:,i);
    end
    
    sigmaSq = zeros(xIdxMax, yIdxMax, s4); % signal strength from variance
    for i = 1 : s4
        for j = 1 : s3
            temp(:,:,j) = filter2(window, imgSeqColor(:, :, j, i).*...
                imgSeqColor(:, :, j, i), 'valid') - lMuSq(:,:,i);
        end
        sigmaSq(:,:,i) = mean(temp, 3);   
    end
    sigma = sqrt( max( sigmaSq, 0 ) );
    ed = sigma * sqrt( wSize^2 * s3 ) + 0.001; % signal strengh
    
% computing weighing map 
    muMap =  exp( -.5 * ( (gMu - .5).^2 /gSig.^2 +  (lMu - .5).^2 /lSig.^2 ) ); % mean intensity weighting map
    normalizer = sum(muMap, 3);
    muMap = muMap ./ repmat(normalizer,[1, 1, s4]);   

    sMap = ed.^p; % signal structure weighting map
    sMap = sMap + 0.001;
    normalizer = sum(sMap,3);
    sMap = sMap ./ repmat(normalizer,[1, 1, s4]);
    
    maxEd = max(ed, [], 3); %  desired signal strength
      
% main loop for motion aware fusion
    fI = zeros(s1, s2, s3); 
    countMap = zeros(s1, s2, s3); 
    countWindow = ones(wSize, wSize, s3);
    xIdx = 1 : stepSize : xIdxMax;
    xIdx = [xIdx xIdx(end)+1 : xIdxMax];
    yIdx = 1 : stepSize : yIdxMax;
    yIdx = [yIdx yIdx(end)+1 : yIdxMax];

    offset = wSize-1;
    for row = 1 : length(xIdx)
        for col = 1 : length(yIdx)
            i = xIdx(row);
            j = yIdx(col);
            blocks = imgSeqColor(i:i+offset, j:j+offset, :, :);
            rBlock = zeros(wSize, wSize, s3);
            for k = 1 : s4
                rBlock = rBlock  + sMap(i, j, k) * ( blocks(:,:,:,k) - lMu(i, j, k) ) / ed(i, j, k);
            end
            if norm(rBlock(:)) > 0
                rBlock = rBlock / norm(rBlock(:)) * maxEd(i, j);
            end
            rBlock = rBlock + sum( muMap(i, j, :) .* lMu(i, j, :) ); 
            fI(i:i+offset, j:j+offset, :) = fI(i:i+offset, j:j+offset, :) + rBlock;
            countMap(i:i+offset, j:j+offset, :) = countMap(i:i+offset, j:j+offset, :) + countWindow;
        end
    end
    fI = fI ./ countMap;
    fI(fI > 1) = 1;
    fI(fI < 0) = 0;
end


