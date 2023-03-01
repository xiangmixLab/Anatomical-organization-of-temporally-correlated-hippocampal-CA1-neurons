

function [place_cells,all_infoScore,all_infoScore_norm] = permutingSpike_adapt_laps_040221(all_neuron,all_behav,temp,occThresh,binsize,all_neuron_num,dir_sign,infosign,coh_thresh)

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

%% null score
[fr_all_dir,ct_all_dir,ctime_all_dir] = calculatinglinearTrackRateMaps_laps(all_neuron,all_behav,binsize,temp,all_neuron_num);
firingrateAll=fr_all_dir{dir_sign}{1}; % 
countAll=ct_all_dir{dir_sign}{1}; % 
countTime=ctime_all_dir{dir_sign}{1}; % 

infoPerSecondnull = zeros(length(firingrateAll),1);
infoPerSpikenull = infoPerSecondnull;


for j = 1:length(firingrateAll)
    MeanFiringRateAll= sum(sum(countAll{j}))/sum(sum(countTime));
    if isempty(firingrateAll{j})
        continue;
    end    
    frt=firingrateAll{j};
%     frt=nansum(frt,1);
%     countTime1=nansum(countTime,1);
    [infoPerSecondnull(j), infoPerSpikenull(j),pxO] = Doug_spatialInfo_020921(frt,countTime,occThresh,MeanFiringRateAll);  
end

%% null coherence

coherencenull = zeros(length(firingrateAll),1);

for j = 1:length(firingrateAll)
    frt=firingrateAll{j};
    [~,coherencenull(j)] = spatial_coherence(frt);  
end

%% shuffled score

infoPerSecondboot = zeros(length(firingrateAll),nboot);
infoPerSpikeboot = infoPerSecondboot;

rng default % for reproduction
trunk_num=10;
for nE = 1:nboot

    all_neuronboot=trunk_shuffle_data_laps(all_neuron,trunk_num);
    [fr_all_dir,ct_all_dir,~] = calculatinglinearTrackRateMaps_laps(all_neuronboot,all_behav,binsize,temp,all_neuron_num);
    firingrateAllt=fr_all_dir{dir_sign}{1}; % 1: left, 2: right
    countAllt=ct_all_dir{dir_sign}{1}; % 1: left, 2: right

    infoPerSecondbootT = zeros(length(firingrateAllt),1);
    infoPerSpikebootT = zeros(length(firingrateAllt),1);

    for j = 1:length(firingrateAllt)
        if isempty(firingrateAllt{j})
            continue;
        end   
        MeanFiringRateAll= sum(sum(countAllt{1,j}))/sum(sum(countTime));
        frt=firingrateAllt{j};
%         frt=nansum(frt,1);
%         countTime1=nansum(countTime,1);
        [infoPerSecondbootT(j), infoPerSpikebootT(j)] = Doug_spatialInfo_020921(frt, countTime,occThresh,MeanFiringRateAll);
    end

    infoPerSecondboot(:,nE) = infoPerSecondbootT; 
    infoPerSpikeboot(:,nE) = infoPerSpikebootT;
end

%% output prepare
infoScoreSecondboot = [infoPerSecondnull,infoPerSecondboot];

infoScoreSpikeboot = [infoPerSpikenull,infoPerSpikeboot];


if isequal(infosign,'sec')
    infoScore = infoScoreSecondboot;    
    infoScoreThresh = quantile(infoScore,0.95,2);
    TinfoPerSecond = table([1:length(infoPerSecondnull)]',infoPerSecondnull,infoScoreThresh,'VariableNames',{'neuron','infoScore','thresh'});
    TinfoPerSecond2 = sortrows(TinfoPerSecond,{'infoScore'},{'descend'});
    place_cells = TinfoPerSecond2.neuron(TinfoPerSecond2.infoScore > TinfoPerSecond2.thresh);
    place_cells=intersect(place_cells,find(coherencenull>coherenceThresh));
end
if isequal(infosign,'spk')
    infoScore = infoScoreSpikeboot;
    infoScoreThresh = quantile(infoScore,0.95,2);
    TinfoPerSecond = table([1:length(infoPerSpikenull)]',infoPerSpikenull,infoScoreThresh,'VariableNames',{'neuron','infoScore','thresh'});
    TinfoPerSecond2 = sortrows(TinfoPerSecond,{'infoScore'},{'descend'});
    place_cells = TinfoPerSecond2.neuron(TinfoPerSecond2.infoScore > TinfoPerSecond2.thresh);
    place_cells=intersect(place_cells,find(coherencenull>coherenceThresh));
end
if isequal(infosign,'all')
    infoScore1 = infoScoreSecondboot;      
    infoScore2 = infoScoreSpikeboot; 
    infoScoreThresh1 = quantile(infoScore1,0.95,2);
    infoScoreThresh2 = quantile(infoScore2,0.95,2);
    
    place_cells = {find(infoPerSecondnull>infoScoreThresh1),find(infoPerSpikenull>infoScoreThresh2)};
    place_cells={intersect(place_cells{1},find(coherencenull>coherenceThresh)),intersect(place_cells{2},find(coherencenull>coherenceThresh))};
end


all_infoScore=[infoPerSecondnull,infoPerSpikenull];
all_infoScore_norm=[(infoPerSecondnull-nanmean(infoScoreSecondboot,2))./nanstd(infoScoreSecondboot,[],2),(infoPerSpikenull-nanmean(infoScoreSpikeboot,2))./nanstd(infoScoreSpikeboot,[],2)];