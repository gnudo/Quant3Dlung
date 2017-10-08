% function for calculating and plotting kernel density estimations
%--------------------------------------------------------------------------
% Date: 2016-07-07
% Author: Goran Lovric
% License: GPL 3 (see LICENSE file in root folder)
%--------------------------------------------------------------------------
clear;clc;close all;
tic

data_files = {...
    '/home/lovric_g/DATA_HD/quant-paper/10x/mE06_mph50_mpp05_mw15/05_ridged_threshed_LocThk.raw'%,...
    %'/home/lovric_g/DATA_HD/quant-paper/10x/mouseE03.csv',...
    %'/home/lovric_g/DATA_HD/quant-paper/10x/mouseE06.csv'
    };

xyz = [256 256 256];

fid = fopen(data_files{1},'r','b');
img = fread(fid, [xyz(1)*xyz(2)*xyz(3)],'float');

img = img(img ~= 0);

[f,xi] = ksdensity(img);
[f1,xi1] = ksdensity(img,'bandwidth',0.1);
[f2,xi2] = ksdensity(img,'bandwidth',0.3);
[f3,xi3] = ksdensity(img,'bandwidth',0.05);

figure;
    plot(xi,f)
    hold on
    plot(xi1,f1)
    plot(xi2,f2)
    plot(xi3,f3)
    hold off
    xlim([0 180])
    xlabel('thicknesses')
    ylabel('density')
    legend('auto bw','bw=0.1','bw=0.3','bw=0.05');
    
toc