%%%%%%% This code is for a single subject data analysis.
%%%%%%% For multiple subjects, you may want to add a for-loop over
%%%%%%% 'subjectid'.
%%%%%%% This code generates FC values for 35 ROI-pairs defined in
%%%%%%% 'ROI_pairs_aal.csv'. Among those 35 pairs, we used only 14 pairs.
%%%%%%% Their pair numbers are 1,8,14,15,16,18,24,25,28,29,30,32,33,35.

%% SubjectID = 1000, 1001, 1002, 1003, 1004, 1005 are healthy controls
%% SubjectID = 2000, 2001, 2002, 2003, 2004, 2005 are MDD

clear all

% Need NIfTI_20140122 tool and spm12
% Also need a few matlabfunctions in matlabfunctions folder included
addpath '/PATH_TO_NIfTI/NIfTI_20140122/';
addpath '/PATH_TO_spm12/spm12/';
addpath '/PATH_TO_matlabfunction/matlabFunctions/';

%% read ROI_pairs_aal.csv file

pairtext = csvread('/PATH_to_ROI_Pairs_AAL_CSV_File/ROI_pairs_aal.csv');
Spairtext = size(pairtext);

% finde the unique aal rois
uniaal = unique(pairtext);

% read all nii file

aal = load_nii('/PATH_to_AAL_nii_FILE/AAL.nii');
aalimg = aal.img;

num_ROI=2;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the subject ID here
subjectid = 1000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data = load_nii(sprintf('/PATH_to_DATA_FILE/%d.nii', subjectid));
nonsmooth_data = data.img;
Sdata = size(nonsmooth_data);
dim = Sdata(1:3);

smooth_data = zeros(Sdata);

% wavelet parameters
waveTname = 'haar';
wavePname = 'rbio3.1';
N.level = 1;
T = Sdata(4);

WT = wavedec3(ones(dim(1), dim(2),dim(3)), N.level ,wavePname);
sWT = size(WT.dec{1});
N.wavecoef = prod(sWT);
coefband = prod(sWT);
usecoefband = 1:coefband;

[TC,TS] = wavedec(zeros(1,T),1,waveTname);
Ntime = length(TC);
waveband = 1:(Ntime) ;

waveT = zeros(N.wavecoef, Ntime, 8);
waveC = zeros(T,coefband, 8);

for t = 1:Sdata(4)
    
    
    temp3d = nonsmooth_data(:,:,:,t);
    smooth3d = imgaussfilt3(temp3d, 1.5); % sigma = 1.5
    smooth_data(:,:,:,t) = smooth3d;
    
    subjbrain = reshape(nonsmooth_data(:,:,:,t), dim(1), dim(2),dim(3));
    
end



count=0;
result = cell( 70 ,4);

for pair=1:Spairtext(1)
    
    % roi number in aal
    aalroi1 = pairtext(pair,1);
    aalroi2 = pairtext(pair,2);
    
    roi1 = (aalimg == aalroi1);
    roi2 = (aalimg == aalroi2);
    
    rclust = cell(num_ROI,1);
    rclust{1} = findn(roi1);
    rclust{2} = findn(roi2);
    
    % get square ROI first
    
    for clust = 1: 2
        temp = rclust{clust};
        
        meancor1 = [round(mean(temp(:,1))) round(mean(temp(:,2))) round(mean(temp(:,3)))];
        totalvoxel = length(temp);
        
        roisize = (1:floor(totalvoxel^(1/3) + 0.001)) -  round(totalvoxel^(1/3)/2+ 0.001) ;
        
        newtemp = zeros( length(roisize)^3 , 3);
        countroi = 1;
        for i=roisize
            for j = roisize
                for k = roisize
                    newtemp(countroi,:) = meancor1 + [ i j k ];
                    countroi= countroi+1;
                end
            end
        end
        
        rclust{clust} = newtemp;
    end
    
    
    % mean voxel analysis
    ROI_timeseries = cell(num_ROI,1);
    
    % get roi time series for mean approach
    for clust = 1:num_ROI
        tempclust = rclust{clust};
        seriesfolder1 = [];
        for i=1:size(tempclust,1)
            temp1 = tempclust(i,:) ;
            tempx1 = temp1(1);
            tempy1 = temp1(2);
            tempz1 = temp1(3);
            tempseri1 = smooth_data(tempx1,tempy1,tempz1,:);
            tempseri1 = tempseri1(:)';
            seriesfolder1 =  [seriesfolder1; tempseri1];
        end
        seriesfolder2 = seriesfolder1;
        % delete any rows (time series at a voxel) contains only zero
        seriesfolder2(~any(seriesfolder1,2),:) = [];
        %tempclust(~any(seriesfolder1,2),:) = [];
        %rclust{clust} = tempclust;
        
        ROI_timeseriesN{clust}  = normseries(seriesfolder2);
    end
    
    
    m1 = mean( ROI_timeseriesN{1},1 );
    m2 = mean( ROI_timeseriesN{2},1 );
    count = count+1;
    result{count,1} = subjectid;
    result{count,2} = 'mean';
    result{count,3} = pair;
    result{count,4} = corr(m1',m2');
    
    
    % double wavelet approach
    
    Sdata = size(nonsmooth_data);
    
    wavedata = zeros(num_ROI, 8,2, T/2);
    wavevar = zeros(num_ROI, 8,2);
    waveform = zeros(num_ROI, 8,2);
    
    for clust = 1:num_ROI
        
        sizeroi = round(length(rclust{clust})^(1/3));
        
        WT = wavedec3(ones(sizeroi, sizeroi,sizeroi), N.level ,wavePname);
        [TC,TS] = wavedec(zeros(1,T),1,waveTname);
        
        N.wavecoef = zeros(8,1);
        
        sWT = size(WT.dec{1});
        N.wavecoef = prod(sWT );
        coefband = prod(sWT);
        usecoefband = 1:coefband;
        
        Ntime = length(TC);
        waveband = 1:(Ntime) ;
        
        
        
        % double wavelet transform
        
        % get roi timeseries
        tempclust = rclust{clust};
        seriesfolder1 = [];
        for i=1:size(tempclust,1)
            temp1 = tempclust(i,:) ;
            tempx1 = temp1(1);
            tempy1 = temp1(2);
            tempz1 = temp1(3);
            tempseri1 = nonsmooth_data(tempx1,tempy1,tempz1,:);
            tempseri1 = tempseri1(:)';
            seriesfolder1 =  [seriesfolder1; tempseri1];
        end
        seriesfolder2 = seriesfolder1;
        % delete any rows (time series at a voxel) contains only zero
        %seriesfolder2(~any(seriesfolder1,2),:) = [];
        %tempclust(~any(seriesfolder1,2),:) = [];
        %rclust{clust} = tempclust;
        
        ROI_time = normseries(double(seriesfolder2));
        ROI_time( find(isnan(ROI_time)) )  = 0;
        
        data = zeros(sizeroi,sizeroi,sizeroi,T);
        
        roicount=1;
        
        for i=1:sizeroi
            for j = 1:sizeroi
                for k = 1:sizeroi
                    data(i,j,k,:) = ROI_time(roicount,:);
                    roicount= roicount+1;
                end
            end
        end
        
        
        waveT = zeros(N.wavecoef, Ntime);
        waveC = zeros(T,coefband);
        
        for loc = 1:8
            
            for t = 1:T
                subjbrain = reshape(data(:,:,:,t), sizeroi, sizeroi,sizeroi);
                WT_temp = wavedec3(subjbrain, N.level ,wavePname);
                waveC(t,:) = reshape(WT_temp.dec{loc},coefband,1 );
            end
            
            for coef = 1:N.wavecoef
                
                subjT = reshape(waveC(:,coef),1,T);
                [waveT(coef, :),TS] = wavedec(subjT,1,waveTname);
                
            end
            
            waveT(~any(waveT,2),:) = [];
            
            t1 = waveT(:,1:(T/2));
            t2 = waveT(:,(1:(T/2))+T/2);
            
            wavedata(clust, loc, 1 ,  :) = mean(t1,1);
            wavedata(clust, loc, 2 ,  :) = mean(t2,1);
            
            ST = size(t1);
            
            vectorT1 = reshape(t1,prod(ST),1  );
            vectorT2 = reshape(t2,prod(ST),1  );
            
            wavevar(clust, loc,1) = var(vectorT1);
            wavevar(clust, loc,2) = var(vectorT2);
            
            waveform(clust, loc,1) = var(vectorT1.* vectorT1);
            waveform(clust, loc,2) = var(vectorT2.* vectorT2);
            
        end
        
    end
    
    chunkcor = zeros(16,1);
    var_weight = zeros(16,1);
    waveform_weight = zeros(16,1);
    
    count_chunk=1;
    
    for spa = 1:8
        for tem = 1:2
            v1 = reshape(wavedata(1,spa,tem,:),T/2,1);
            v2 = reshape(wavedata(2,spa,tem,:),T/2,1);
            
            chunkcor(count_chunk) = corr(v1,v2);
            
            var_weight(count_chunk) = mean(wavevar(:, spa,tem));
            waveform_weight(count_chunk) = mean(waveform(:, spa,tem));
            count_chunk = count_chunk+1;
        end
    end
    
    cor_waveform = chunkcor .* waveform_weight / sum(waveform_weight);
    cor_variance = chunkcor .* var_weight / sum(var_weight);
    
    % 'dw_waveform' is the functional connectivity based on the double
    % wavelet approach
    count = count+1;
    result{count,1} = subjectid;
    result{count,2} = 'dw_waveform';
    result{count,3} = pair;
    result{count,4} = sum(cor_waveform);
    
    % 'dw_variance' will not be used in our analysis
    count = count+1;
    result{count,1} = subjectid;
    result{count,2} = 'dw_variance';
    result{count,3} = pair;
    result{count,4} = sum(cor_variance);
    
    
end


cd('PATH_to_FOLDER_to_SAVE_RESULTS/FOLDER_NAME')
cell2csv(sprintf('%d.csv', subjectid) ,result)
