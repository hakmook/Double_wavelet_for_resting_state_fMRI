
% Multiple subjest resting state data
% Double-wavelet approach for 'gau100' datatype, i.e., all the voxels have
% the same resting state fMRI data

% In this sample code, the underyling true correlation is 0.


% Need NIfTI_20140122 tool and spm12
% Also need a few matlabfunctions in matlabfunctions folder included

addpath '/PATH_TO_NIfTI/NIfTI_20140122/';
addpath '/PATH_TO_spm12/spm12/';
addpath '/PATH_TO_matlabfunction/matlabFunctions/';


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




%%%%% For wavelet only

waveTname = 'haar';
wavePname = 'rbio3.1';

N.level=1;

% get 2d wavelet coefficient size
[C, S]  = wavedec2(zeros(dim(1), dim(2)), N.level ,wavePname);

% number of low-low wavelet coefficient
plot.size = S(1,1)*S(1,2);
%plot.size = sum(S(:,1) .* S(:,2));
coefband = 1:length(C);

% number of 2d wavelet coefficient
N.wavecoef = length(C);

% create wavelet coefficient matrix
waveC = zeros( T, length(C));
% 1D wavelet coefficient
[TC,TS] = wavedec(zeros(1,T),1,waveTname);

Ntime = length(TC);
waveband = 1:Ntime ;

waveT = zeros(length(C), length(TC) );

wave_pval_glm = zeros(R,Nsimu);

wavecon = zeros( Nsimu, R,Nsubj, length(coefband));
wavecon_sub = zeros(Nsubj, R);


wavedata =  zeros(R, length(coefband), length(waveband));


result = cell(100, 6);
count = 0;

for cortype = 1 : length(cordata)
    
    for sim = 1: Nsimu
        
        %for i = 1:length(allfile)
        cd('PATH_TO_DATA') % cd to the folder where the data are located
        
        load(sprintf('Exp2_Sim%d_Gau100' ,sim));
        
        for sub = 1:Nsubj
            
            for roi = 1:R
                for t = 1:T
                    subjbrain = reshape(data(sub,roi,:,:,t), dim(1), dim(2));
                    [waveC(t,:), S] = wavedec2(subjbrain, N.level ,wavePname);
                end
                
                for coef = 1:N.wavecoef
                    
                    subjT = reshape(waveC(:,coef),1,T);
                    [waveT(coef, :),TS] = wavedec(subjT,1,waveTname);
                    
                end
                wavedata( roi, :,:) = waveT(coefband,waveband);
            end
            
            sdata = size(wavedata);
            
            spadata1 = (1:(sdata(2)*1/4)) + 0* (sdata(2)*1/4) ;
            spadata2 = (1:(sdata(2)*1/4)) + 1* (sdata(2)*1/4) ;
            spadata3 = (1:(sdata(2)*1/4)) + 2* (sdata(2)*1/4) ;
            spadata4 = (1:(sdata(2)*1/4)) + 3* (sdata(2)*1/4) ;
            
            spadataAll = [spadata1;spadata2;spadata3;spadata4];
            
            tempdataA = (1:(sdata(3)/2)) + 0*(sdata(3)/2);
            tempdataB = (1:(sdata(3)/2)) + 1*(sdata(3)/2);
            
            tempdataAll = [tempdataA; tempdataB];
            
            
            for spa = 1:4
                
                
                for temp = 1:2
                    
                    spadata = spadataAll(spa,:);
                    tempdata = tempdataAll(temp,:);
                    
                    d1 = reshape( mean( wavedata(1,spadata,tempdata),2),length(tempdata),1);
                    d2 = reshape( mean( wavedata(2,spadata,tempdata),2),length(tempdata),1);
                    
                    wavecor(spa, temp) = corr(d1, d2);
                    
                    var1 = var(reshape(wavedata(1,spadata,tempdata), length(spadata)*length(tempdata)  ,1));
                    var2 = var(reshape(wavedata(2,spadata,tempdata), length(spadata)*length(tempdata)  ,1));
                    
                    variance(spa, temp) = mean([var1 var2])  ;
                    
                    
                    vector1 = reshape(wavedata(1,spadata,tempdata), length(spadata)*length(tempdata)  ,1);
                    vector2 = reshape(wavedata(2,spadata,tempdata), length(spadata)*length(tempdata)  ,1);
                    
                    variance_waveform(spa, temp) = mean([var(vector1.*vector1) var(vector2.*vector2)])  ;
                    
                end
                
            end
            
            tempvar = sum( reshape(variance(:,:), 1,8) );
            compcor = sum(sum(reshape(variance(:, :), 4,2) .* reshape(wavecor(:, :), 4,2)))/tempvar;
            
            tempvar_waveform = sum( reshape(variance_waveform(:,:), 1,8) );
            compcor_waveform = sum(sum(reshape(variance_waveform(:, :), 4,2) .* reshape(wavecor(:, :), 4,2)))/tempvar_waveform;
            
            count=count+1;
            result{count,1} = 'gau100';
            result{count,2} = sim;
            result{count,3} = '0';
            result{count,4} = sub;
            result{count,5} = compcor;
            result{count,6} = 'dw_variance';
            
            count=count+1;
            result{count,1} = 'gau100';
            result{count,2} = sim;
            result{count,3} = '0';
            result{count,4} = sub;
            result{count,5} = compcor_waveform;
            result{count,6} = 'dw_waveform';
            
            
        end
        
        cell2csv( sprintf('dw_gau100_%d.csv', sim) ,result)
        
    end
end








