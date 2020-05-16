function result = MSEF(images_rgb)
[H, W, ~, N] = size(images_rgb);
%% calculate adaptive well-exposedness weight maps, namesly weight_maps_1
% compute luminance
images_luminance = zeros(H, W, N); % initialize
for i = 1 : N
    images_ycbcr = rgb2ycbcr(images_rgb(:, :, :, i));
    images_luminance(:, :, i) = images_ycbcr(:, :, 1);
    %figure,imshow(images_rgb(:, :, :, i));
end
% compute mean value of luminance
mean_value = mean(mean(images_luminance));
mean_value = reshape(mean_value, N, 1);
% set the sigma value
sigma = 0.2; % 0.2
% compute adaptive well-exposedness weight maps denoted as weight_maps_1
weight_maps_1 = zeros(H, W, N); % initialize
for j = 1 : N
    weight_maps_1(:, :, j) = exp(-0.5 * (images_luminance(:, :, j) - (1 - mean_value(j))).^2 /sigma /sigma);
end
% calculate 3-D color gradient, namely weight_maps_2
weight_maps_2 = zeros(H, W, N);
for k = 1 : N
    [weight_maps_2(:, :, k), ~, ~] = color_gradient(images_rgb(:, :, :, k));
end
% initial weight weight_maps
w1 = 1.0; w2 = 2.2;
weight_maps = (weight_maps_1.^w1) .* (weight_maps_2.^w2);
% refine weight maps
weight_refine = zeros(H, W, N);% initialize
for m = 1 : N
    weight_refine(:, :, m) = imgaussfilt(weight_maps(:, :, m), 3);
end
% normalize weight maps
weight_refine = weight_refine + 1e-12; %avoids division by zero
weight_refine = weight_refine./repmat(sum(weight_refine, 3),[1 1 N]);
% determine the number of decomposition layers
if N > 3
    nlev = 7;
else
    nlev = 8;
end
% create empty pyramid
pyr = gaussian_pyramid(zeros(H, W, 3), nlev);
%% multi-scale exposure fusion
for i = 1 : N
    % construct pyramid
	pyrW = gaussian_pyramid(weight_refine(:,:,i), nlev);
	pyrI = laplacian_pyramid(images_rgb(:,:,:,i), nlev);
    % carry out fusion
    for l = 1 : nlev
        w = repmat(pyrW{l}, [1, 1, 3]);
        pyr{l} = pyr{l} + w.*pyrI{l};
    end
end
%% reconstruct
result = reconstruct_laplacian_pyramid(pyr);
%% pyramid fusion
result = uint8(255*result);
end

    
