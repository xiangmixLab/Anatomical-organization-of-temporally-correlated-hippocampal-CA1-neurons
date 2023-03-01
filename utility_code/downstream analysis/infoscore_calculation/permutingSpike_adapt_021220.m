
function [place_cells,all_infoScore,all_infoScore_norm] = permutingSpike_adapt_021220(neuron,behavpos,behavtime,temp,occThresh,binsize,small_velo,infosign)

if ~exist('occThresh','var') || isempty(occThresh)
%     occThresh = 0.2; %061219 adjust code
      occThresh = 0.1;
end

if ~exist('temp','var') || isempty(temp)
    temp = 'S';
end

if ~exist('binsize','var') || isempty(binsize)
    binsize=15;
end

% trim low fr neurons, they tend to have high infoscore but that doesn't
% necessarily mean they are place cells...

%% prepare statistics
nboot = 100;
countTimeThresh = [0.1 inf]; 



maxbehavROI=[0 0 max(behavpos(:,1)),max(behavpos(:,2))];

neuron.S=C_to_peakS(neuron.C);
thresh=0.1*max(neuron.C,[],2);

[behav_struct.position,behav_struct.time]=align_neuron_behav(neuron,behavtime,behavpos); % resampled position for MI
velo_vec=behav_velo_cal(behav_struct,[],'r');

%% null score
[firingrateAll,countAll,~,countTime,~,~,bininfo,cpt] = calculatingCellSpatialForSingleData_Suoqin(neuron,behavpos,behavtime,maxbehavROI,binsize,1:size(neuron.C,1),thresh,temp,[],[],countTimeThresh,small_velo);

infoPerSecondnull = zeros(length(firingrateAll),1);
infoPerSpikenull = infoPerSecondnull;
MInull=infoPerSecondnull;

for j = 1:length(firingrateAll)
    MeanFiringRateAll= sum(sum(countAll{j}))/sum(sum(countTime));
    if isempty(firingrateAll{j})
        continue;
    end    
    frt=firingrateAll{j};
    [infoPerSecondnull(j), infoPerSpikenull(j),pxO] = Doug_spatialInfo_020921(frt,cpt,occThresh,MeanFiringRateAll);  
    
    nSt=neuron.S(j,:);
    nSt(nSt<0.1*nanmax(nSt))=0;
    nSt=nSt>0;
%     [MInull(j), posterior, occupancy_map, prob_being_active, likelihood] = extract_2D_information(nSt, behav_struct.position, bininfo.xpos, bininfo.ypos, velo_vec>small_velo);
    nSt1=nSt;
    nSt1(velo_vec<=small_velo)=0;
    tic;
    [MInull(j)] = extract_2D_information_adapted(nSt1, countAll{j},cpt);
    toc;
end

%% shuffled score
S0 = neuron.S;
C0 = neuron.C;
time0 = neuron.time;

infoPerSecondboot = zeros(length(firingrateAll),nboot);
infoPerSpikeboot = infoPerSecondboot;
MIboot = infoPerSecondboot;

rng default % for reproduction
% rotate_leng=floor(size(S0,2)/(nboot+1));
% C_dat=trunk_split_data(C0,100,'col');
% S_dat=trunk_split_data(S0,100,'col');
firingrateAllt_collection={};
deltaT=[0.01:(0.99-0.01)/nboot:0.99];
for nE = 1:nboot
%     neuronboot=neuron.copy;
%     C1=trunk_shuffle_data_pre_split(C_dat);
%     S1=trunk_shuffle_data_pre_split(S_dat);
%     neuronboot.S = S1;
%     neuronboot.C = C1;
%     neuronboot=trunk_shuffle_data_neuron_simple(C_dat,S_dat,time0);
%     neuronboot=trunk_rotate_data_neuron_simple(C0,S0,time0,nE*(size(S0,2)/15-100)/nboot);
    rotate_leng=floor(deltaT(nE)*size(C0,2));
    S1=trunk_rotate_data(S0,rotate_leng,'right');
    C1=trunk_rotate_data(C0,rotate_leng,'right');
    neuronboot.S = S1;
    neuronboot.C = C1;
    neuronboot.time = time0;
    [firingrateAllt,countAllt] = calculatingCellSpatialForSingleData_Suoqin(neuronboot,behavpos,behavtime,maxbehavROI,binsize,1:size(neuronboot.C,1),thresh,temp,[],[],countTimeThresh,small_velo);

    infoPerSecondbootT = zeros(length(firingrateAllt),1);
    infoPerSpikebootT = zeros(length(firingrateAllt),1);
    MIbootT = zeros(length(firingrateAllt),1);

    for j = 1:length(firingrateAllt)
        if isempty(firingrateAll{j})
            continue;
        end   
        MeanFiringRateAll= sum(sum(countAllt{1,j}))/sum(sum(countTime));
        frt=firingrateAllt{j};
%         frt=filter2DMatrices(firingrateAllt{j},1);
%         frt=imgaussfilt(firingrateAll{j},1);
        [infoPerSecondbootT(j), infoPerSpikebootT(j)] = Doug_spatialInfo_020921(frt, countTime,occThresh,MeanFiringRateAll);
        
        nSt=neuronboot.S(j,:);
        nSt(nSt<0.1*nanmax(nSt))=0;
%         [MIbootT(j), posterior, occupancy_map, prob_being_active, likelihood] = extract_2D_information(nSt>0, behav_struct.position, bininfo.xpos, bininfo.ypos, velo_vec>small_velo);
        nSt=nSt>0;
        nSt1=nSt;
        nSt1(velo_vec<=small_velo)=0;
        [MIbootT(j)] = extract_2D_information_adapted(nSt1, countAll{j},cpt);
    end

    infoPerSecondboot(:,nE) = infoPerSecondbootT; 
    infoPerSpikeboot(:,nE) = infoPerSpikebootT;
    MIboot(:,nE) = MIbootT;
    firingrateAllt_collection{nE}=firingrateAllt;
end

%% output prepare
infoScoreSecondboot = [infoPerSecondnull,infoPerSecondboot];

infoScoreSpikeboot = [infoPerSpikenull,infoPerSpikeboot];

MIboot = [MInull,MIboot];


if isequal(infosign,'sec')
    infoScore = infoScoreSecondboot;    
%     infoScore(neuron_lowFR,:) = [];
%     infoScoreThresh = quantile(infoScore(:),0.95);
    infoScoreThresh = quantile(infoScore,0.95,2);
%     infoScoreThresh = quantile(infoScore,0.99,2);% change to 0.99 061219
    TinfoPerSecond = table([1:length(infoPerSecondnull)]',infoPerSecondnull,infoScoreThresh,'VariableNames',{'neuron','infoScore','thresh'});
%     TinfoPerSecond(neuron_lowFR,:) = [];
    TinfoPerSecond2 = sortrows(TinfoPerSecond,{'infoScore'},{'descend'});
    place_cells = TinfoPerSecond2.neuron(TinfoPerSecond2.infoScore > TinfoPerSecond2.thresh);
end
if isequal(infosign,'spk')
    infoScore = infoScoreSpikeboot;
%     infoScore(neuron_lowFR,:) = [];
%     infoScoreThresh = quantile(infoScore(:),0.95);
    infoScoreThresh = quantile(infoScore,0.95,2);
%     infoScoreThresh = quantile(infoScore,0.99,2);% change to 0.99 061219
    TinfoPerSecond = table([1:length(infoPerSpikenull)]',infoPerSpikenull,infoScoreThresh,'VariableNames',{'neuron','infoScore','thresh'});
%     TinfoPerSecond(neuron_lowFR,:) = [];
    TinfoPerSecond2 = sortrows(TinfoPerSecond,{'infoScore'},{'descend'});
    place_cells = TinfoPerSecond2.neuron(TinfoPerSecond2.infoScore > TinfoPerSecond2.thresh);
end
if isequal(infosign,'all')
    infoScore1 = infoScoreSecondboot;      
    infoScore2 = infoScoreSpikeboot; 
%     infoScore(neuron_lowFR,:) = [];
%     infoScoreThresh = quantile(infoScore(:),0.95);
    infoScoreThresh1 = quantile(infoScore1,0.95,2);
    infoScoreThresh2 = quantile(infoScore2,0.95,2);
    infoScoreThresh3 = quantile(MIboot,0.95,2);
%     infoScoreThresh = quantile(infoScore,0.99,2);% change to 0.99 061219
    
    place_cells = {find(infoPerSecondnull>infoScoreThresh1),find(infoPerSpikenull>infoScoreThresh2),find(MInull>infoScoreThresh3)};
end


all_infoScore=[infoPerSecondnull,infoPerSpikenull,MInull];
all_infoScore_norm=[(infoPerSecondnull-nanmean(infoScoreSecondboot,2))./nanstd(infoScoreSecondboot,[],2),(infoPerSpikenull-nanmean(infoScoreSpikeboot,2))./nanstd(infoScoreSpikeboot,[],2),(MInull-nanmean(MIboot,2))./nanstd(MIboot,[],2)];