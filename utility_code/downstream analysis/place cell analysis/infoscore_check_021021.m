function [all_infoscore] = infoscore_check_021021(neuron,behavpos,behavtime,temp,occThresh,binsize,small_velo)

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

countTimeThresh = [0.1 inf]; 

maxbehavROI=[0 0 max(behavpos(:,1)),max(behavpos(:,2))];

neuron.S=C_to_peakS(neuron.C);
thresh=0.1*max(neuron.C,[],2);

[firingrateAll,countAll,~,countTime,~,~,bininfo,cpt] = calculatingCellSpatialForSingleData_040321(neuron,behavpos,behavtime,maxbehavROI,binsize,1:size(neuron.C,1),thresh,temp,[],[],countTimeThresh,small_velo);

infoPerSecondnull = zeros(length(firingrateAll),1);
infoPerSpikenull = infoPerSecondnull;

for j = 1:length(firingrateAll)
    MeanFiringRateAll= sum(sum(countAll{j}))/sum(sum(countTime));
    if isempty(firingrateAll{j})
        continue;
    end
    
    frt=firingrateAll{j};
    [infoPerSecondnull(j), infoPerSpikenull(j)] = Doug_spatialInfo_020921(frt,countTime,occThresh,MeanFiringRateAll);
end
all_infoscore=[infoPerSecondnull,infoPerSpikenull];
% nboot=100;
% S0 = neuron.S;
% C0 = neuron.C;
% time0 = neuron.time;
% 
% infoPerSecondboot = zeros(length(firingrateAll),nboot);
% infoPerSpikeboot = infoPerSecondboot;
% 
% rng default % for reproduction
% C_dat=trunk_split_data(C0,100,'col');
% S_dat=trunk_split_data(S0,100,'col');
% for nE = 1:nboot
%     neuronboot=trunk_shuffle_data_neuron_simple(C_dat,S_dat,time0);
% 
%     [firingrateAllt,countAllt,~,countTimet] = calculatingCellSpatialForSingleData_Suoqin(neuronboot,behavpos,behavtime,maxbehavROI,binsize,1:size(neuronboot.C,1),thresh,temp,[],[],countTimeThresh,small_velo);
% 
%     infoPerSecondbootT = zeros(length(firingrateAllt),1);infoPerSpikebootT = zeros(length(firingrateAllt),1);
% 
%     for j = 1:length(firingrateAllt)
%         MeanFiringRateAll= sum(sum(countAllt{1,j}))/sum(sum(countTimet));
%         [infoPerSecondbootT(j), infoPerSpikebootT(j)] = Doug_spatialInfo_020921(frt, countTimet,occThresh,MeanFiringRateAll);
%     end
% 
%     infoPerSecondboot(:,nE) = infoPerSecondbootT; 
%     infoPerSpikeboot(:,nE) = infoPerSpikebootT;
% end
% 
% all_infoscore=[infoPerSecondnull,infoPerSpikenull];
% all_infoscore_norm=[infoPerSecondnull_norm,infoPerSpikenull_norm];
% 
% infoScoreSecondboot = [infoPerSecondnull,infoPerSecondboot];
% infoScoreSpikeboot = [infoPerSpikenull,infoPerSpikeboot];
% 
% infoPerSecondnull=(infoPerSecondnull-mean(infoScoreSecondboot,2))./std(infoScoreSecondboot,[],2); % z-score 
% infoPerSpikenull=(infoPerSpikenull-mean(infoScoreSpikeboot,2))./std(infoScoreSpikeboot,[],2); % z-score 
