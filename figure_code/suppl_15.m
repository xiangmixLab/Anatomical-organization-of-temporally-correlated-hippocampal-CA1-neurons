% SUPPL 18: lap-by-lap stability for linear track experiment

foldername=foldername_LT;
%% 1. lap-by-lap calculation
for i=1:length(foldername)
    neuronName=[foldername{i},'\','neuronIndividuals_new.mat'];
    behavName={
        [foldername{i},'\','Horizontal1_Behav.mat'];
        [foldername{i},'\','Horizontal2_Behav.mat'];
        [foldername{i},'\','Vertical_Behav.mat'];
    };
    load(neuronName);
    for j=1:length(neuronIndividuals_new)
        load(behavName{j});
        [nnew{i,j},behavnew{i,j},laps_dir{i,j}]=lineartrack_laps_separate(neuronIndividuals_new{j},behav);
        
    end
end
[nnew,behavnew]=rearrange_laps(nnew,behavnew);
% [nnew,behavnew]=laps_avg_neuroDat_040622(nnew,behavnew); % AVG direction, neuron and behav relationship

%% 2. lap firing field
[fr_all,ctime_all] = calculatinglinearTrackRateMaps_laps_noDir(nnew,behavnew,10,'S');


%% 3. lap-lap corr
meanCorrOdd={};
meanCorrEven={};
meanCorr={};
lap_ratemap={}
for i=1:6
    for j=1:3
        
        cellNum=length(fr_all{i,j}{1,1});
        if cellNum==0
            cellNum=length(fr_all{i,j}{2,1});
        end
        
        for k=1:cellNum
            lap_ratemap{i,j}{k,1}=track_lap_ratemap_by_cell(fr_all{i,j},k);
            
            corrVal1=corrcoef(lap_ratemap{i,j}{k,1}{1}');
%             for p=1:size(corrVal1,1) % previous codes also not involve
%             diagonal deletion
%                 corrVal1(p,p)= nan;
%             end
            
            meanCorrOdd{i,j}(k,1)=mean(corrVal1(:),'omitnan');
            
            corrVal2=corrcoef(lap_ratemap{i,j}{k,1}{2}');
%             for p=1:size(corrVal2,1)
%                 corrVal2(p,p)= nan;
%             end
            
            meanCorrEven{i,j}(k,1)=mean(corrVal2(:),'omitnan');
        end
        meanCorr{i,j}=max([meanCorrOdd{i,j},meanCorrEven{i,j}],[],2);
    end
end



%% lap-lap corr, pc-z
% calculate pc
all_behav={};
for i=1:length(foldername)
    load([foldername{i},'\','behav.mat']);
    for j=1:size(behavIndividuals,2) 
        all_behav{i,j}=behavIndividuals{j};
        all_behav{i,j}.VidObj=[];
    end
end

all_pc_z=cell(6,3);
all_infoscore_z=cell(6,3);
all_infoscore_norm_z=cell(6,3);
all_coherence_z=cell(6,3);
tic;
for i=1:6
    load([foldername{i},'\','neuronIndividuals_new.mat']);
    for j=1:3
    
        [place_cells,infoScore,infoScore_norm,coherencee] = permutingSpike_adapt_041222_LT(neuronIndividuals_new{j},all_behav{i,j}.position,all_behav{i,j}.time,'S',0,10,5,'all',0.4);  

        all_pc_z{i,j}=place_cells;
        all_infoscore_z{i,j}=infoScore;
        all_infoscore_norm_z{i,j}=infoScore_norm;
        all_coherence_z{i,j}=coherencee;

    end
end

% pc-z ratemaps
meanCorrOdd_PCz={};
meanCorrEven_PCz={};
meanCorr_PCz={};

for i=1:6
    for j=1:3
        
        cellIdx=all_pc_z{i,j}{2};
        
        for k=1:length(cellIdx)
            lap_ratemap{i,j}{k,1}=track_lap_ratemap_by_cell(fr_all{i,j},cellIdx(k));
            
            corrVal1=corrcoef(lap_ratemap{i,j}{k,1}{1}');

            meanCorrOdd_PCz{i,j}(k,1)=mean(corrVal1(:),'omitnan');
            
            corrVal2=corrcoef(lap_ratemap{i,j}{k,1}{2}');
            meanCorrEven_PCz{i,j}(k,1)=mean(corrVal2(:),'omitnan');
        end
        meanCorr_PCz{i,j}=max([meanCorrOdd_PCz{i,j},meanCorrEven_PCz{i,j}],[],2);
    end
end


%% llustration
% ratemaps, hor1 high/low

%%%%%%%%%%%%%%%%%%%%%% test ratemaps
% odd
for mice=1:6
    lmaps=lap_ratemap{mice,1};
    figure;
    ctt=1;
    for k=1:length(lmaps)
         if (meanCorrOdd{mice,1}(k)>=0.6)&&(meanCorrOdd{mice,1}(k)<=0.9)
             subplot(8,8,ctt)
             imagesc(lmaps{k}{1});
             title(num2str(k))
             ctt=ctt+1;
         end
    end
end
% even
for mice=1:6
    lmaps=lap_ratemap{mice,1};
    figure;
    ctt=1;
    for k=1:length(lmaps)
         if (meanCorrEven{mice,1}(k)>=0.6)&&(meanCorrEven{mice,1}(k)<=0.9)
             subplot(8,8,ctt)
             imagesc(lmaps{k}{2});
             title(num2str(k))
             ctt=ctt+1;
         end
    end
end
%%%%%%%%%%%%%%%%%%%%%% test ratemaps

%A
lmaps=lap_ratemap{4,1};
figure;
subplot(221)
imagesc(lmaps{53}{1});colorbar
title([num2str(meanCorrOdd{4,1}(53)),',',num2str(max(lmaps{53}{1}(:)))])
subplot(222)
imagesc(lmaps{150}{1});colorbar
title([num2str(meanCorrOdd{4,1}(150)),',',num2str(max(lmaps{150}{1}(:)))])

subplot(223)
imagesc(lmaps{136}{1});colorbar
title([num2str(meanCorrOdd{4,1}(136)),',',num2str(max(lmaps{136}{1}(:)))])
subplot(224)
imagesc(lmaps{228}{2});colorbar
title([num2str(meanCorrEven{4,1}(228)),',',num2str(max(lmaps{228}{2}(:)))])


%C
figure;
lmaps=lap_ratemap{5,2};
subplot(221)
imagesc(lmaps{98}{1});colorbar
title([num2str(meanCorrOdd{5,2}(98)),',',num2str(max(lmaps{98}{1}(:)))])
subplot(222)
imagesc(lmaps{209}{1});colorbar
title([num2str(meanCorrOdd{5,2}(209)),',',num2str(max(lmaps{209}{1}(:)))])

subplot(223)
imagesc(lmaps{93}{1});colorbar
title([num2str(meanCorrOdd{5,2}(93)),',',num2str(max(lmaps{93}{1}(:)))])
subplot(224)
imagesc(lmaps{165}{1});colorbar
title([num2str(meanCorrOdd{5,2}(165)),',',num2str(max(lmaps{165}{1}(:)))])


%E
figure;
lmaps=lap_ratemap{2,3};

subplot(221)
imagesc(lmaps{272}{2});colorbar
title([num2str(meanCorrEven{2,3}(272)),',',num2str(max(lmaps{272}{2}(:)))])
subplot(222)
imagesc(lmaps{312}{2});colorbar
title([num2str(meanCorrEven{2,3}(312)),',',num2str(max(lmaps{312}{2}(:)))])

subplot(223)
imagesc(lmaps{41}{2});colorbar
title([num2str(meanCorrEven{2,3}(41)),',',num2str(max(lmaps{41}{2}(:)))])
subplot(224)
imagesc(lmaps{274}{2});colorbar
title([num2str(meanCorrEven{2,3}(274)),',',num2str(max(lmaps{274}{2}(:)))])



% B,D,F lap-lap corr, allcells/ pcs, hor1
subplot(131)
c1=cell2mat(meanCorr(:,1));
c2=cell2mat(meanCorr_PCz(:,1));
cdfplot(c1);hold on;
cdfplot(c2);
nanmedian(c1)
nanmedian(c2)


subplot(132)
c3=cell2mat(meanCorr(:,2));
c4=cell2mat(meanCorr_PCz(:,2));
cdfplot(c3);hold on;
cdfplot(c4);
nanmedian(c3)
nanmedian(c4)

subplot(133)
c5=cell2mat(meanCorr(:,3));
c6=cell2mat(meanCorr_PCz(:,3));
cdfplot(c5);hold on;
cdfplot(c6);
nanmedian(c5)
nanmedian(c6)

c1=distribution_based_resampling(c1,10,10);
c2=distribution_based_resampling(c2,10,10);
c3=distribution_based_resampling(c3,10,10);
c4=distribution_based_resampling(c4,10,10);
c5=distribution_based_resampling(c5,10,10);
c6=distribution_based_resampling(c6,10,10);

% pc vs all
% [h,p1]=kstest2([mean(meanCorr{1,1});mean(meanCorr{2,1});mean(meanCorr{3,1});mean(meanCorr{4,1});mean(meanCorr{5,1});mean(meanCorr{6,1})],[mean(meanCorr_PCz{1,1});mean(meanCorr_PCz{2,1});mean(meanCorr_PCz{3,1});mean(meanCorr_PCz{4,1});mean(meanCorr_PCz{5,1});mean(meanCorr_PCz{6,1})]);
% [h,p2]=kstest2([mean(meanCorr{1,2});mean(meanCorr{2,2});mean(meanCorr{3,2});mean(meanCorr{4,2});mean(meanCorr{5,2});mean(meanCorr{6,2})],[mean(meanCorr_PCz{1,2});mean(meanCorr_PCz{2,2});mean(meanCorr_PCz{3,2});mean(meanCorr_PCz{4,2});mean(meanCorr_PCz{5,2});mean(meanCorr_PCz{6,2})]);
% [h,p3]=kstest2([mean(meanCorr{1,3});mean(meanCorr{2,3});mean(meanCorr{3,3});mean(meanCorr{4,3});mean(meanCorr{5,3});mean(meanCorr{6,3})],[mean(meanCorr_PCz{1,3});mean(meanCorr_PCz{2,3});mean(meanCorr_PCz{3,3});mean(meanCorr_PCz{4,3});mean(meanCorr_PCz{5,3});mean(meanCorr_PCz{6,3})]);
[h,p1]=kstest2(c1,c2);
[h,p2]=kstest2(c3,c4);
[h,p3]=kstest2(c5,c6);

%% first half and second half corr
corrRes={};
for i=1:length(foldername)
    neuronName=[foldername{i},'\','neuronIndividuals_new.mat'];
    behavName={
        [foldername{i},'\','Horizontal1_Behav.mat'];
        [foldername{i},'\','Horizontal2_Behav.mat'];
        [foldername{i},'\','Vertical_Behav.mat'];
    };

    load(neuronName);
    all_behav={};
    for j=1:length(behavName)
        load(behavName{j});
        all_behav{j}=behav;
    end
    corrRes(i,:)=LT_firstSecondHalf_corr_dir(neuronIndividuals_new, all_behav, 10)';
end

corrRes_mat_odd=[]
corrRes_mat_even=[]
for i=1:6
    for j=1:3
        corrRes_mat_odd(i,j)=nanmean(corrRes{i,j}(:,1));
        corrRes_mat_even(i,j)=nanmean(corrRes{i,j}(:,2));
    end
end