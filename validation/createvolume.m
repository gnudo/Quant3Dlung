% 
%--------------------------------------------------------------------------
% Date: 2017-07-12
% Author: Goran Lovric
%--------------------------------------------------------------------------
clear;clc;close all;


%--------------------------------------------------------------------------
% 1.) Constants
%--------------------------------------------------------------------------

resolution = 512;       % [px] number of pixels in each 3D direction
n_balls = 350;          % [1]  total number of balls within the volume
dataname = 'fullValidation';
R_balls  = [10 15 20 25 30 35 50];       % [px] different ball radii
p_balls  = [5 10 15 25 30 10 5];         % [%] percentages of ball counts
                                         % in distribution

%--------------------------------------------------------------------------
% 2.) Random distribution
%--------------------------------------------------------------------------
B_mat = zeros(resolution,resolution,resolution);
kk = 1

while kk < n_balls+1;
    MC_fac = rand(1);
    
    R_ind = 1+floor( rand(1)*length(R_balls) );
    R = R_balls(R_ind);
    p_1 = p_balls(R_ind)/100;
    
    if MC_fac > p_1
        continue
    end
    
    r_pos = round (1 + rand(1,3) .* (resolution-1) );
    
    quit = false;
    
    R_mod = R + 4;   % for increasing distance between balls and borders
    
    % (1) We test whether the ball is within the volume
    if r_pos(1) - R_mod < 1 || r_pos(2) - R_mod < 1 || r_pos(3) - R_mod < 1
        continue
    end
	
	if r_pos(1) + R_mod > resolution || r_pos(2) + R_mod > resolution || r_pos(3) + R_mod > resolution
        continue
    end
    
    % (2) We test whether we have space to place the ball in the matrix    
    for ii = -R_mod:R_mod
        for jj = -R_mod:R_mod
            for ll = -R_mod:R_mod
                xx = r_pos(1)+ii;
                yy = r_pos(2)+jj;
                zz = r_pos(3)+ll;
                if sqrt (ii.^2 + jj.^2 + ll.^2) > R_mod;
                    continue
                end
                if B_mat(xx,yy,zz) == 1;
                    quit = true;
                    disp('Overlapping balls!')
                    break;
                end
            end
            if quit
                break;
            end
        end
        if quit
            break;
        end
    end
    
    % (3) If we have space, then we place the balls it
    if quit == false;
        for ii = -R:R%(R-1)
            for jj = -R:R%(R-1)
                for ll = -R:R%(R-1)
                xx = r_pos(1)+ii;
                yy = r_pos(2)+jj;
                zz = r_pos(3)+ll;
                if sqrt (ii.^2 + jj.^2 + ll.^2) < R;
                    B_mat(xx,yy,zz) = 1;
                end
                end
            end
        end
        R_vec(kk) = R;
        kk = kk + 1
    end
end


%--------------------------------------------------------------------------
% 3.) Folder manipulation
%--------------------------------------------------------------------------
target_dir = dataname;

if isdir(target_dir)
    rmdir(target_dir,'s'); % delete previously created files
end
if exist(target_dir,'dir') ~= 7
    mkdir(target_dir);
end


%--------------------------------------------------------------------------
% 4.) Statistical analysis
%--------------------------------------------------------------------------
xvalues = min(R_vec):max(R_vec);
[nelements,centers] = hist(R_vec,xvalues);
bar(centers,nelements)

A = [floor(centers); nelements];

fileID = fopen('quant_data.txt','w');
fprintf(fileID,'%6s %12s\n','R [px]','Amount');
fprintf(fileID,'%6.2f %12.0f\n',A);
fclose(fileID);


%--------------------------------------------------------------------------
% 5.) Save images
%--------------------------------------------------------------------------
zspace = '0000';
for ind = 1:size(B_mat,3)
    zer=zspace(1:end-length(num2str(ind)));
	file=sprintf('%s/%s%d.tif',target_dir,zer,ind)
	img = uint8( (2^8-1).*(B_mat(:,:,ind)) ); % save as 8bit
	imwrite(img,file,'tif');
end
