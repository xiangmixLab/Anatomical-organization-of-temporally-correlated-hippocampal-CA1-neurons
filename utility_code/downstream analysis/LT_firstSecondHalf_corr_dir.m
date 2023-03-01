function corrRes=LT_firstSecondHalf_corr_dir(neuronIndividuals_new, all_behav, binsize)

nnew={};
behavnew={};
laps_dir={};
for j=1:length(neuronIndividuals_new)
    [nnew{1,j},behavnew{1,j},laps_dir{1,j}]=lineartrack_laps_separate(neuronIndividuals_new{j},all_behav{j});
end

[nnew,behavnew]=rearrange_laps(nnew,behavnew);
% [nnew,behavnew]=laps_avg_neuroDat_040622(nnew,behavnew); % AVG direction, neuron and behav relationship

%% 2. lap firing field   

[fr_all,ctime_all] = calculatinglinearTrackRateMaps_laps_noDir(nnew,behavnew,binsize,'S');

%% 3. first half and second half

corrRes={};
for i=1:length(neuronIndividuals_new)
    frt=fr_all{i};
    frt_dir1=frt(:,1);
    frt_dir2=frt(:,2);
    
    % clear empty cells
    frt_dir1=frt_dir1(~cellfun('isempty',frt_dir1));
    frt_dir2=frt_dir2(~cellfun('isempty',frt_dir2));
    
    % half1
    frt_dir1_h1=frt_dir1(1:floor(size(frt_dir1,1)/2));
    frt_dir1_h2=frt_dir1(floor(size(frt_dir1,1)/2)+1:end);

    
    frt_dir1_h1_cellMat=LT_half_arrange(frt_dir1_h1);
    frt_dir1_h2_cellMat=LT_half_arrange(frt_dir1_h2);
        
    % half2
    frt_dir2_h1=frt_dir1(1:floor(size(frt_dir2,1)/2));
    frt_dir2_h2=frt_dir1(floor(size(frt_dir2,1)/2)+1:end);
        
    frt_dir2_h1_cellMat=LT_half_arrange(frt_dir2_h1);
    frt_dir2_h2_cellMat=LT_half_arrange(frt_dir2_h2);
    
    % corr
    frt_dir1_corr=[];
    frt_dir2_corr=[];
    for k=1:length(frt_dir1_h1_cellMat)
        try
            frt1=corrcoef(frt_dir1_h1_cellMat{k},frt_dir1_h2_cellMat{k});
            frt_dir1_corr(k,1)=frt1(2);
        catch
            frt_dir1_corr(k,1)=nan;
        end
        try
            frt2=corrcoef(frt_dir2_h1_cellMat{k},frt_dir2_h2_cellMat{k});
            frt_dir2_corr(k,1)=frt2(2);
        catch
            frt_dir2_corr(k,1)=nan; 
        end
    end
    
    corrRes{i,1}=[frt_dir1_corr,frt_dir2_corr];
end
