function connected_component(datasegment)
 
files = dir(fullfile(datasegment,'*.tif'));
 
img = imread(fullfile(datasegment,files(1).name));
 
y = zeros(size(img,1),size(img,2),length(files),'uint8');
y(:,:,1) = img;
 
parfor k = 2:length(files)
    img = imread(fullfile(datasegment,files(k).name));
    y(:,:,k) = img;
end
 
CC = bwconncomp(y,6);
L = labelmatrix(CC);
 
a = histcounts(L(:),unique(L(:)));
[~,a2] = sort(a);
 
tissue_id = a2(end-1)-1;
img_seg = zeros(size(y),'uint8');
img_seg(L==tissue_id) = 1;
 
write_dir = strcat(datasegment,'/labeled');
mkdir(write_dir)
 
parfor k = 1:length(files)
    
imwrite(logical(img_seg(:,:,k)),fullfile(write_dir,files(k).name),'tif','Compression','none');
end