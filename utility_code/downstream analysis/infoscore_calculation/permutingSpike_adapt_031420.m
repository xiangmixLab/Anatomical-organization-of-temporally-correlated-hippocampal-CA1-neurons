
function [place_cells,all_infoScore,all_infoScore_norm] = permutingSpike_adapt_031420(neuron,behavpos,behavtime,temp,occThresh,binsize,small_velo,infosign)

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

%% null score
[firingrateAll,countAll,~,countTime,~,~,bininfo,cpt] = calculatingCellSpatialForSingleData_040321(neuron,behavpos,behavtime,maxbehavROI,binsize,1:size(neuron.C,1),thresh,temp,[],[],countTimeThresh,small_velo);

infoPerSecondnull = zeros(length(firingrateAll),1);
infoPerSpikenull = infoPerSecondnull;

for j = 1:length(firingrateAll)
    MeanFiringRateAll= sum(sum(countAll{j}))/sum(sum(countTime));
    if isempty(firingrateAll{j})
        continue;
    end    
    frt=firingrateAll{j};
    [infoPerSecondnull(j), infoPerSpikenull(j),pxO] = Doug_spatialInfo_020921(frt,cpt,occThresh,MeanFiringRateAll);  
end

%% shuffled score
S0 = neuron.S;
C0 = neuron.C;
time0 = neuron.time;

infoPerSecondboot = zeros(length(firingrateAll),nboot);
infoPerSpikeboot = infoPerSecondboot;

rng default % for reproduction
firingrateAllt_collection={};
deltaT=[0.01:(0.99-0.01)/nboot:0.99];
for nE = 1:nboot
    rotate_leng=floor(deltaT(nE)*size(C0,2));
    S1=trunk_rotate_data(S0,rotate_leng,'right');
    C1=trunk_rotate_data(C0,rotate_leng,'right');
    neuronboot.S = S1;
    neuronboot.C = C1;
    neuronboot.time = time0;
    [firingrateAllt,countAllt] = calculatingCellSpatialForSingleData_040321(neuronboot,behavpos,behavtime,maxbehavROI,binsize,1:size(neuronboot.C,1),thresh,temp,[],[],countTimeThresh,small_velo);

    infoPerSecondbootT = zeros(length(firingrateAllt),1);
    infoPerSpikebootT = zeros(length(firingrateAllt),1);

    for j = 1:length(firingrateAllt)
        if isempty(firingrateAll{j})
            continue;
        end   
        MeanFiringRateAll= sum(sum(countAllt{1,j}))/sum(sum(countTime));
        frt=firingrateAllt{j};
        [infoPerSecondbootT(j), infoPerSpikebootT(j)] = Doug_spatialInfo_020921(frt, countTime,occThresh,MeanFiringRateAll);        
    end

    infoPerSecondboot(:,nE) = infoPerSecondbootT; 
    infoPerSpikeboot(:,nE) = infoPerSpikebootT;
    firingrateAllt_collection{nE}=firingrateAllt;
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
end
if isequal(infosign,'spk')
    infoScore = infoScoreSpikeboot;
    infoScoreThresh = quantile(infoScore,0.95,2);
    TinfoPerSecond = table([1:length(infoPerSpikenull)]',infoPerSpikenull,infoScoreThresh,'VariableNames',{'neuron','infoScore','thresh'});
    TinfoPerSecond2 = sortrows(TinfoPerSecond,{'infoScore'},{'descend'});
    place_cells = TinfoPerSecond2.neuron(TinfoPerSecond2.infoScore > TinfoPerSecond2.thresh);
end
if isequal(infosign,'all')
    infoScore1 = infoScoreSecondboot;      
    infoScore2 = infoScoreSpikeboot; 
    infoScoreThresh1 = quantile(infoScore1,0.95,2);
    infoScoreThresh2 = quantile(infoScore2,0.95,2);
    place_cells = {find(infoPerSecondnull>infoScoreThresh1),find(infoPerSpikenull>infoScoreThresh2)};

end


all_infoScore=[infoPerSecondnull,infoPerSpikenull];
all_infoScore_norm=[(infoPerSecondnull-nanmean(infoScoreSecondboot,2))./nanstd(infoScoreSecondboot,[],2),(infoPerSpikenull-nanmean(infoScoreSpikeboot,2))./nanstd(infoScoreSpikeboot,[],2)];