clear;

%wenjianlujing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
path = './images/';
files = dir(path);

length_file = length(files);
a = files.name;
for i = 3:length_file
    fprintf('current set is %f\n',i/10*10);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
temp_path = [path,files(i).name];

II = load_images(temp_path);
%%%resort
[h,w,c,n]= size(II);
for ii = 1:n
    meann(ii) = mean(mean(mean(II(:,:,:,ii))));
end
[mean2, index] = sort(meann);

for jj=1:n
    I(:,:,:,jj) = II(:,:,:,index(jj));
end

dst(:,:,1,:) = 0.299*I(:,:,1,:) + 0.587*I(:,:,2,:) + 0.114*I(:,:,3,:);
dst(:,:,2,:) = -0.147*I(:,:,1,:) - 0.289*I(:,:,2,:) + 0.436*I(:,:,3,:);
dst(:,:,3,:) = 0.615*I(:,:,1,:) - 0.515*I(:,:,2,:)- 0.100*I(:,:,3,:);
[height,width,channel]=size(I(:,:,:,1));
clear temp;

tic
R=code_for_CSVT20191(dst,[1 1 1]);
runtime = toc;
FR(:,:,1)= R(:,:,1) + 1.14*R(:,:,3);  
FR(:,:,2)= R(:,:,1) - 0.39*R(:,:,2) - 0.58*R(:,:,3);   
FR(:,:,3)= R(:,:,1) + 2.03*R(:,:,2);

fprintf('the running time is %f\n', runtime);
imshow(FR);
% write_path_root = '/media/king/AC88794B887914D4/TIP_1011/20181020/code_test/';
% write_path = [write_path_root,files(i).name,'.png'];
% 
% imwrite(FR,write_path);
clear height;
clear width;
clear I;
clear dst;
clear FR;
clear R;
clear II;
clear n;
clear index;
clear meann;
end