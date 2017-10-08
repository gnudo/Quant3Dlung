function segmenter(datasource,minpeakheight,minpeakprominence,maxwidth,v)
warning ('off','all');

datadestin = [datasource(1:end-1), '_mph', minpeakheight, '_mpp', ...
              minpeakprominence, '_mw', maxwidth, '/'];

files = dir([datasource '*.tif']);

minpeakheight = eval(minpeakheight);
minpeakprominence = eval(minpeakprominence);
maxwidth = eval(maxwidth);

if ~exist('v')
    v = 'false';
end    

struct_element_bg = 7;%20;
if ~isdir(datadestin)
    mkdir(datadestin); % delete previously created files
end
if strcmp(lower(v),'true')
    datadestin_01 = [datadestin '01_bg_corrected/'];
    datadestin_02 = [datadestin '02_original_otsu/'];
    datadestin_03 = [datadestin '03_smaller_otsu/'];
    datadestin_04 = [datadestin '04_ridged_image/'];
    datadestin_05 = [datadestin '05_ridged_threshed/'];
    mkdir(datadestin_01);
    mkdir(datadestin_02);
    mkdir(datadestin_03);
    mkdir(datadestin_04);
    mkdir(datadestin_05);
else
    datadestin_01 = datadestin;  % because of parfor
    datadestin_02 = datadestin;  % because of parfor
    datadestin_03 = datadestin;  % because of parfor
    datadestin_04 = datadestin;  % because of parfor
    datadestin_05 = datadestin;
end

tic
parfor n=1:length(files)
zspace = '0000';
zer=zspace(1:end-length(num2str(n)));

%--------------------------------------------------------------------------
% (0) Load data
%--------------------------------------------------------------------------
img_path = [datasource files(n).name];
img_orig1 = imread(img_path);

%--------------------------------------------------------------------------
% (1) Equalize background
%--------------------------------------------------------------------------

% % (a) CLAHE
% img_newww = least_square_plane(double(img_orig1));
% img_newww = adapthisteq(img_orig1,'NumTiles',[64 64]);

% (b) Gaussian background subtraction
% h = fspecial('gaussian',100,0.5);
% img_newww = imfilter(img_orig1,h,'symmetric','conv');
% img_orig1 = double(img_orig1)-1.2.*double(img_newww)+mean(img_orig1(:));

% (c) morphological background subtraction and plane fitting
mmm = size(img_orig1,1); 
nnn = size(img_orig1,2); 

[YY, XX] = meshgrid(1:nnn, 1:mmm);
z1 = double(imerode(img_orig1,strel('disk',struct_element_bg)));
[muhat,sigmahat] = normfit(z1(:));
z1(z1<muhat-3*sigmahat | z1>muhat+3*sigmahat) = NaN;

X = [ones(mmm*nnn, 1), YY(:), XX(:)];
X1 = X(~isnan(z1(:)),:);

M = X1\z1(~isnan(z1(:)));
back = reshape(X*M, mmm, nnn);
img_orig1 = double(img_orig1)-back;
img_orig1 = img_orig1 - min(img_orig1(:));
img_orig1 = uint8(img_orig1);

if strcmp(lower(v),'true')
    file_bg = sprintf('%s%s%d.tif',datadestin_01,zer,n);
    imwrite(img_orig1,file_bg,'tif','Compression','none');
end

%--------------------------------------------------------------------------
% (2) Transform into a quadratic image
%--------------------------------------------------------------------------
img_orig = zeros(max(size(img_orig1)),'uint8');
img_orig(1:size(img_orig1,1),1:size(img_orig1,2)) = img_orig1;

size_img_orig = size(img_orig);
img_ridged = img_orig;

%--------------------------------------------------------------------------
% (3) Otsu Thresholding of original image
%--------------------------------------------------------------------------
thresh_level = graythresh(uint8(img_orig1));
img_orig_seg = im2bw(img_orig(1:size(img_orig1,1),1:size(img_orig1,2)),thresh_level);

if strcmp(lower(v),'true')
    file_bg = sprintf('%s%s%d.tif',datadestin_02,zer,n);
    imwrite(logical(img_orig_seg),file_bg,'tif','Compression','none');
    
    % Save the reduced Otsu
    thresh_level_red = thresh_level*0.75;
    img_orig_seg_red = im2bw(img_orig(1:size(img_orig1,1),1:size(img_orig1,2)),thresh_level_red);
    file_bg = sprintf('%s%s%d.tif',datadestin_03,zer,n);
    imwrite(logical(img_orig_seg_red),file_bg,'tif','Compression','none');
    
end

%--------------------------------------------------------------------------
% (4) Ridge enhencer in all 4x directions
%--------------------------------------------------------------------------
% (1) Horizontal direction
for ii=1:size(img_orig,1)
    profile = img_orig(ii,:);
    [pks,locs,w,p] = findpeaks(double(profile),'MinPeakHeight', minpeakheight, ...
                               'MinPeakProminence',minpeakprominence,'Annotate','extents');
    locs = locs(w<maxwidth);
    img_ridged(ii,locs) = 255;
end

% % (2) Vertical direction
for jj=1:size(img_orig,2)
    profile = img_orig(:,jj);
    [pks,locs,w,p] = findpeaks(double(profile),'MinPeakHeight', minpeakheight, ...
                              'MinPeakProminence',minpeakprominence,'Annotate','extents');
%     locs = locs(w<maxwidth);
    img_ridged(locs,jj) = 255;
end

% (3) 45 degree 
for jj=-(size(img_orig,1)-3):(size(img_orig,2)-3)
    profile = diag(img_orig,jj);
    [pks,locs,w,p] = findpeaks(double(profile),'MinPeakHeight', minpeakheight, ...
                              'MinPeakProminence',minpeakprominence,'Annotate','extents');
	locs = locs(w<maxwidth);
    profile(locs) = 255;
	img1 = diag(profile,jj);
    img_ridged(img1==255)=255;
end

% (4) 135 degree 
img_orig = rot90(img_orig);
img_ridged = rot90(img_ridged);
for jj=-(size(img_orig,1)-3):(size(img_orig,2)-3)
    profile = diag(img_orig,jj);
    [pks,locs,w,p] = findpeaks(double(profile),'MinPeakHeight', minpeakheight, ...
                               'MinPeakProminence',minpeakprominence,'Annotate','extents');
    locs = locs(w<maxwidth);
    profile(locs) = 255;
	img1 = diag(profile,jj);
    img_ridged(img1==255)=255;
end
img_orig=rot90(img_orig,3);
img_ridged=rot90(img_ridged,3);

img_ridged = img_ridged(1:size(img_orig1,1),1:size(img_orig1,2));

if strcmp(lower(v),'true')
    file_bg = sprintf('%s%s%d.tif',datadestin_04,zer,n);
    imwrite(img_ridged,file_bg,'tif','Compression','none');
end
%--------------------------------------------------------------------------
% (5) Otsu Thresholding and noise removal
%--------------------------------------------------------------------------
img_seg = zeros(size(img_orig1),'uint8');
img_seg(img_ridged > thresh_level*255) = 1;

%--------------------------------------------------------------------------
% (6) Morphological operations
%--------------------------------------------------------------------------
img_seg = morpho_operation(double(img_seg),30);

%--------------------------------------------------------------------------
% (7) Write image file
%--------------------------------------------------------------------------
file=sprintf('%s%s%d.tif',datadestin_05,zer,n);
imwrite(logical(img_seg),file,'tif','Compression','none');
end
toc