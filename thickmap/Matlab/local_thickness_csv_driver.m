% function for creating CSV-histogram files from the thickness map. The CSV
% file is created at the same location as the RAW file
%--------------------------------------------------------------------------
% Date: 2016-07-05
% Author: Yannis Vogiatzis with modifications from Goran Lovric
% License: GPL 3 (see LICENSE file in root folder)
%--------------------------------------------------------------------------
function local_thickness_csv_driver(raw_file,x,y,z)

x = eval(x);
y = eval(y);
z = eval(z);

% raw_file = '/home/lovric_g/DATA_HD/quant-paper/10x/mE02_mph50_mpp05_mw15/05_ridged_threshed_mini2_LocThk.raw';
% xyz = [64 64 64];
xyz = [x y z];


diameters = [0:200];

% directory = directories{l};

saveFileName = strcat(raw_file(1:end-4),'.mat');
cor_coeff = 0.07968;

fid = fopen(raw_file,'r','b');
img = fread(fid, [xyz(1)*xyz(2)*xyz(3)],'float');

thickness_map = reshape(img,[xyz(1) xyz(2) xyz(3)]);
thickness_map = permute(thickness_map,[2 1 3]);

diameter_counts = hist(thickness_map(:),diameters);

% neglecting the background value
%diameters = diameters(2:end);

vol_struct = struct('diameter_counts',diameter_counts,'diameters',diameters);    

parsave(saveFileName);
csvwrite(strcat(raw_file(1:end-4),'.csv'),[vol_struct.diameters.' , vol_struct.diameter_counts.'])

end