% function y_details = quad_decom_Fast_WLS( vec_h, vec_v, beta, lambda, T)
function Zinit = quad_decom_Fast_WLS(Zinit, VFx,VFy, beta, LAMBDA, T)

[H, W] = size(Zinit);
lambda = zeros(T,1);
for t=1:T
    lambda(t) = 1.5* 4^(T-t)/(4^T-1)*LAMBDA;
end
epsilon = 2;%0.0001; % epsilon is set to 2
%%%parameters along x direction 
row = zeros(1,W);
a = zeros(1,W);
b = zeros(1,W);
c = zeros(1,W);
d = zeros(1,W);
lambda_tx = zeros(1,W);
%%%parameters along y direction
col = zeros(1,H);
aa = zeros(1,H);
bb = zeros(1,H);
cc = zeros(1,H);
dd = zeros(1,H);
lambda_ty = zeros(1,H);
    
for t=1:T %iteration for separable optimization
%%%definition of matrices coefficients  

    for h = 1: H % for each row, A is a matrix of WxW
        for i=1:W-1
            lambda_tx(i) = lambda(t)/(abs(VFx(h,i))^beta + epsilon);
        end
        a(1) = -lambda_tx(1);
        b(1) = 1+lambda_tx(1); 
        c(1) = a(1);
        d(1) = Zinit(h, 1)-lambda_tx(1)* VFx(h,1);        
        for i = 2:W-1
            a(i) = -lambda_tx(i);
            b(i) = 1+lambda_tx(i)+lambda_tx(i-1);
            c(i) = a(i);
            d(i) = Zinit(h, i) +lambda_tx(i-1)*VFx(h,i-1) - lambda_tx(i)*VFx(h,i);          
        end
        b(W) = 1+lambda_tx(W-1);
        d(W) = Zinit(h, W) +lambda_tx(W-1)* VFx(h,W-1);     
       %%%%1D solver along x direction                        
        row(1) = d(1)/b(1);
                     
    %%%%setting values to '\tilde b'
        for i = 2:W
            b(i) = b(i) - a(i-1)*c(i-1)/b(i-1);
            row(i) = (d(i) - a(i-1)*row(i-1))/b(i);
        end
            
        for i=(W-1):-1:1
            row(i) = row(i) - c(i)*row(i+1)/b(i);
        end
             
        Zinit(h,:) = row(:);
    end
%%%Definition of matrix coefficients
      
    for w = 1:W % for each col, A is a matrix is HxH
        for i=1:H-1
            lambda_ty(i) = lambda(t)/(abs(VFy(i,w))^beta + epsilon);
        end
        aa(1) = - lambda_ty(1);
        bb(1) = 1+ lambda_ty(1);
        cc(1) = aa(1);
        dd(1) = Zinit(1, w) - lambda_ty(1)* VFy(1,w);
        for i = 2:H-1
            aa(i) = -lambda_ty(i);
            bb(i) = 1+lambda_ty(i)+lambda_ty(i-1);
            cc(i) = aa(i);
            dd(i) = Zinit(i, w) + lambda_ty(i-1)*VFy(i-1,w) - lambda_ty(i)*VFy(i,w); 
        end        
        bb(H) = 1+ lambda_ty(H-1); 
        dd(H) = Zinit(H, w) + lambda_ty(H-1)* VFy(H-1,w);
  %%%1D solver along y direction
                        
        col(1) = dd(1)/bb(1);
               
            %setting values to '\tilde b'
        for i = 2:H
            bb(i) = bb(i) - aa(i-1)*cc(i-1)/bb(i-1);
            col(i) = (dd(i) - aa(i-1)*col(i-1))/bb(i);
        end
            
        for i=(H-1):-1:1
            col(i) = col(i) - cc(i)*col(i+1)/bb(i);
        end
            
        Zinit(:,w) = col(:);
    end
        
end
    
    
end






