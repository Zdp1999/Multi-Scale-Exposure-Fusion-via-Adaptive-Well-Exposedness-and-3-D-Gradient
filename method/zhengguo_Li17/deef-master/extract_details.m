%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [L_d] = extract_details(I, R_int)
height = size(I, 1);
width = size(I, 2);
ldrNum = size(I, 4);
LOG2_LUT = zeros(256,1); 

for ii=1:256
    LOG2_LUT(ii) = log2(ii);
end

T_W = zeros(256,1);
for ii=1:128
    T_W(ii) = ii; %%1;
end
for ii=129:256
    T_W(ii) = 257-ii; %%1;
end

%%%Convert the input images into log domain
Y_int = zeros(height,width);
Y_ldr = zeros(height,width,ldrNum);
Y_tmp = zeros(height,width,ldrNum);
for ii=1:height
    for jj=1:width
        for k=1:ldrNum
            Y_tmp(ii,jj,k) = floor((77*I(ii,jj,1,k)+150*I(ii,jj,2,k)+29*I(ii,jj,3,k)+128)/256);
            Y_ldr(ii,jj,k) = LOG2_LUT(Y_tmp(ii,jj,k)+1);
        end
        Y_int(ii,jj) = LOG2_LUT(floor((77*R_int(ii,jj,1)+150*R_int(ii,jj,2)+29*R_int(ii,jj,3)+384)/256));
    end
end

%%%generation of vector field via spetral fusion
gvf_h = zeros(height, width);
gvf_v = zeros(height, width);
gvf_Y_h = zeros(height, width);
gvf_Y_v = zeros(height, width);
Zh = zeros(ldrNum,2);
WWW = zeros(ldrNum,2);
for ii=1:height-1
    for jj=1:width-1
         for k=1:ldrNum
            WWW(k,1) = T_W(Y_tmp(ii,jj,k)+1)*T_W(Y_tmp(ii,jj+1,k)+1);
            WWW(k,2) = T_W(Y_tmp(ii,jj,k)+1)*T_W(Y_tmp(ii+1,jj,k)+1);
        end
        WWW1 = max(WWW(:,1));
        WWW2 = max(WWW(:,2));       
        for k=1:ldrNum
            Zh(k,1) = (Y_ldr(ii,jj+1,k)-Y_ldr(ii,jj,k))*WWW(k,1)/WWW1;
            Zh(k,2) = (Y_ldr(ii+1,jj,k)-Y_ldr(ii,jj,k))*WWW(k,2)/WWW2;
        end
        [Uh, Sh, Vh] = svd(Zh, 0);
        gvf_Y_h(ii,jj) = Y_int(ii,jj+1)-Y_int(ii,jj);
        gvf_Y_v(ii,jj) = Y_int(ii+1,jj)-Y_int(ii,jj);       
        Z_tmp = [gvf_Y_h(ii,jj),gvf_Y_v(ii,jj)]*Sh*Vh./(abs(gvf_Y_h(ii,jj))+abs(gvf_Y_v(ii,jj))+0.000000001);
        gvf_h(ii,jj) = Z_tmp(1); 
        gvf_v(ii,jj) = Z_tmp(2);
    end
     for k=1:ldrNum
        WWW(k,2) = T_W(Y_tmp(ii,width,k)+1)*T_W(Y_tmp(ii+1,width,k)+1);
    end
    WWW2 = max(WWW(:,2));    
    for k=1:ldrNum
        Zh(k,1) = 0;
        Zh(k,2) = (Y_ldr(ii+1,width,k)-Y_ldr(ii,width,k))*WWW(k,2)/WWW2;
    end
    [Uh, Sh, Vh] = svd(Zh, 0);
    gvf_Y_h(ii,width) = 0;
    gvf_Y_v(ii,width) = Y_int(ii+1,width)-Y_int(ii,width);
    Z_tmp = [gvf_Y_h(ii,width),gvf_Y_v(ii,width)]*Sh*Vh./(abs(gvf_Y_v(ii,width))+0.000000001);
    gvf_h(ii,width) = Z_tmp(1);
    gvf_v(ii,width) = Z_tmp(2);
end
for jj=1:width-1
         for k=1:ldrNum
            WWW(k,1) = T_W(Y_tmp(height,jj,k)+1)*T_W(Y_tmp(height,jj+1,k)+1);
        end
        WWW1 = max(WWW(:,1));
    for k=1:ldrNum
        Zh(k,1) = (Y_ldr(height,jj+1,k)-Y_ldr(height,jj,k))*WWW(k,1)/WWW1;
        Zh(k,2) = 0;
    end
    [Uh, Sh, Vh] = svd(Zh, 0);
    gvf_Y_h(height,jj) = Y_int(height,jj+1)-Y_int(height,jj);
    gvf_Y_v(height,jj) = 0;  
    Z_tmp = [gvf_Y_h(height,jj),gvf_Y_v(height,jj)]*Sh*Vh./(abs(gvf_Y_h(height,jj))+0.000000001);
    gvf_h(height,jj) = Z_tmp(1);
    gvf_v(height,jj) = Z_tmp(2);
end

%%%correction of the vector field
 for ii=1:height
    for jj=1:width
        if gvf_h(ii,jj)*gvf_Y_h(ii,jj)<0
            gvf_h(ii,jj) = -gvf_h(ii,jj);
        end
        if gvf_v(ii,jj)*gvf_Y_v(ii,jj)<0
            gvf_v(ii,jj) = -gvf_v(ii,jj);
        end
%         if (gvf_h(ii,jj)*gvf_Y_h(ii,jj)+gvf_v(ii,jj)*gvf_Y_v(ii,jj))<0
%             gvf_h(ii,jj) = -gvf_h(ii,jj);
%             gvf_v(ii,jj) = -gvf_v(ii,jj);
%         end        
   end
end

% extract detals from the vector field
beta = 0.75;%1.25; %0.75;
lambda = 1;  
T = 3;
Zinit = zeros(height,width);
[L_d] = quad_decom_Fast_WLS(Zinit,gvf_h, gvf_v, beta, lambda, T); % fast WLS

% [L_d] = quad_decom_Fast(gvf_h, gvf_v, height, width, beta, lambda); % WLS

% figure('Name', 'details_layer'); imshow(2.^L_d); 
imwrite(2.^L_d, 'fusion_details.png','png');
% figure('Name', 'log_details_layer'); imshow(L_d);
% imwrite(L_d, 'log_details_layer.png','png');
end

