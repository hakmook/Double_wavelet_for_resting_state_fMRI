%%%%%%% This code is for simulate subject data only
%%%%%%%

%% Number of ROIs = 2, T = 256

clear all

% Need NIfTI_20140122 tool and spm12
% Also need a few matlabfunctions in matlabfunctions folder included

addpath '/PATH_TO_NIfTI/NIfTI_20140122/';
addpath '/PATH_TO_spm12/spm12/';
addpath '/PATH_TO_matlabfunction/matlabFunctions/';


Nsubj = 20; % Number of subjects
Nsimu = 500; % number of simulation


% AR(1) coefficient
phi = 0.6;

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
spat_phi = 2;
sigma_spat=1.5;
dist_matrix1 = cell(num_ROI,1);
dx1 = cell(num_ROI,1);
for clust = 1:num_ROI
    dx1{clust} = ones(num_voxROI(clust),1);
    distV = pdist(coord_ROI{clust},'euclidean');
    
    dist_matrix1{clust} = squareform(distV);
end

spat_COV_true = sigma_spat*exp(- (dist_matrix1{1}./spat_phi) );
chol_spat = chol(spat_COV_true);


for sim = 1: Nsimu

    for sub = 1: Nsubj

        beta_temp_spat1 = mvnrnd(zeros(Voxel,1),eye(Voxel));

        beta_spat1 = chol_spat'*beta_temp_spat1';
        
        error = zeros(R, dim(1) ,dim(2), T);

        count_beta = 0;

        for j = 1:dim(2)
            for i = 1:dim(1)
                count_beta = count_beta + 1;

                temp_v1 = mvnrnd(zeros(R,1), SIgma);
               

                for roi = 1:R
                    error(roi,i,j,1) = temp_v1(roi)+ beta_spat1(count_beta);
                    
                end

                for h = 2:T
                    w = mvnrnd(zeros(R,1), SIgma);
                    for roi = 1:R 
                        error(roi,i,j,h) = phi*error(roi,i,j,h-1) + w(roi)+ beta_spat1(count_beta); 
                    end
                end


            end
        end
    
        
    Y_data = zeros(R, 10,10, 256);
    
    for j = 1:R
        test = error(j,:,:,:);
        test1 = normseries(reshape(test, 100, 256)); % normalize series
        Y_data(j,:,:,:) = reshape(test1, 10,10, 256);
    end

      data(sub,:, :,:,:) = Y_data;
 
    end
    
	save(sprintf('Exp2_Sim%d_cortype' ,sim)    ,'data')    
    
    
end





