
% Multiple subjest resting state data
% mean (AVG_FC) approach for 'gau100' datatype, i.e., all the voxels have
% the same resting state fMRI data

% In this sample code, the underyling true correlation is 0.


% Need NIfTI_20140122 tool and spm12
% Also need a few matlabfunctions in matlabfunctions folder included

addpath '/PATH_TO_NIfTI/NIfTI_20140122/';
addpath '/PATH_TO_spm12/spm12/';
addpath '/PATH_TO_matlabfunction/matlabFunctions/';

% read data

T = 256; t = T;
R = 2; num_ROI = 2;
alpha = 0.1;
Voxel = 100;

Nsimu = 500;
Nsubj = 50;
dim = [10 10];


%%%%%%%%%%%%%%%%%%%

corname = {'0'}; % the underlying true correaltion used for data generation
cordata = 0; % the underlying true correaltion used for data generation
datatpye = {'gau100' }; % generated data type in {'exp', 'gau0', 'gau100'}

dname = datatpye;


result = cell(100, 6);
count = 0;

for cortype = 1 : length(cordata)
    
    for sim = 1: Nsimu
        
       cd('PATH_TO_DATA') % cd to the folder where the data are located
        
        load(sprintf('Exp2_Sim%d_Gau100' ,sim));
        
        for sub = 1: Nsubj
            
            smooth_error = zeros(R, dim(1),dim(2),T);
            for roi = 1:R
                for h = 1:T
                    tempY = data(sub,roi,:,:,h);%error(roi,:,:,h);
                    tempY = reshape(tempY,dim(1),dim(2) );
                    smooth_error(roi,:,:,h) = imfilter(tempY,H,'replicate');
                    
                end
            end
            smooth_error = reshape(smooth_error,2, dim(1)*dim(2),T);
            d1 = reshape( mean( smooth_error(1,:,:),2),T,1);
            d2 = reshape( mean( smooth_error(2,:,:),2),T,1);
            meancor = corr(d1, d2);
            
            
            
            count=count+1;
            result{count,1} = 'gau100';
            result{count,2} = sim;
            result{count,3} = '0';
            result{count,4} = sub;
            result{count,5} = meancor;
            result{count,6} = 'mean';
            
        end
        
        
    end
end


cd('PATH_TO_FOLDER_TO_SAVE_RESULTS')
cell2csv( sprintf('mean_gau100_%d.csv', Nsimu), result)



