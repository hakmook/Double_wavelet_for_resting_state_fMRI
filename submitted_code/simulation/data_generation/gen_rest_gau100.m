%%%%%%% This code is for simulate subject data only
%%%%%%%

%% Number of ROIs = 2, T = 256

clear all


% Need NIfTI_20140122 tool and spm12
% Also need a few matlabfunctions in matlabfunctions folder included

addpath '/PATH_TO_NIfTI/NIfTI_20140122/';
addpath '/PATH_TO_spm12/spm12/';
addpath '/PATH_TO_matlabfunction/matlabFunctions/';

Nsubj = 50; % Number of subjects
Nsimu = 500; % number of simulation


% AR(1) coefficient
phi = 0.6;
% phi = 0.4:0.05:0.85;


T = 256; t = 256;
R = 2; num_ROI = 2;
Voxel = 100;

num_run = 1;
num_images = 256; % images per run
num_voxROI = Voxel*ones(num_ROI,1);

% testcor: define the correlation of interest,
% choose one from {0, 0.1, 0.2,..., 1.0}
testcor = 0;

SIgma = [1 , testcor;
    testcor, 1];
rng(testcor*100)


dim = [10, 10];

% Create coordinates for each voxel
coord_ROI = cell(R,1);
temp = [];
for j = 1:dim(1)
    for k = 1:dim(2)
        temp = [temp; j, k];
    end
end
coord_ROI{1} = temp;
coord_ROI{2} = [temp(:,1),temp(:,2)+20];

sqrt_sigma = sqrt(5)*eye(R); % variance at each ROI

data = zeros(Nsubj,R, dim(1), dim(2) , T);


% spatial distance matrix

%%%%%%%%%%%%%%%%%%%
filter_sigma = 100;
%filter_sigma = 1.1:0.1:2;

filtersize = 10;
H = fspecial('gaussian',filtersize, filter_sigma);


sqrt_sigma = sqrt(5)*eye(R); % variance at each ROI





for sim = 1: Nsimu
    
    data = zeros(Nsubj,R, dim(1), dim(2) , T);
    
     for sub = 1: Nsubj
        
        rho_vec = zeros(1, 10);
        error = zeros(R, dim(1), dim(2),T);
        
        
        for j = 1:dim(2)
            for i = 1:dim(1)
                temp_v = mvnrnd(zeros(R,1), SIgma);
                for roi = 1:R
                    error(roi,i,j,1) = temp_v(roi) ;
                end
                
                for h = 2:T
                    w = mvnrnd(zeros(R,1), SIgma);
                    for roi = 1:R
                        error(roi,i,j,h) = phi*error(roi,i,j,h-1) + w(roi) ;
                    end
                end
                
                
            end
        end
        
        
        Y_smoothing = error;
        
        % Create Y's
        
        
        Y_data = zeros(R, 10,10, T);
        for j = 1:R
            test = Y_smoothing(j,:,:,:);
            test1 = normseries(reshape(test, 100, T)); % normalize series
            Y_data(j,:,:,:) = reshape(test1, 10,10, T);
        end
        
        
        
        data(sub,:, :,:,:) = Y_data;
        
    end
    
    save(sprintf('Exp2_Sim%d_Gau100' ,sim), 'data')
    
    
end

