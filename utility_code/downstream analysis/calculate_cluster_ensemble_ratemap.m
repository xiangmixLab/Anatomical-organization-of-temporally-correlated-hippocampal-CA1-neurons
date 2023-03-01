function calculate_cluster_ensemble_ratemap(neuron,idx)

neuron.C=neuron.C(idx,:);
neuron.S=neuron.S(idx,:);

[firingrateAll,countAll,~,countTime,~,~,bininfo,cpt] = calculatingCellSpatialForSingleData_040321(neuron,behavpos,behavtime,maxbehavROI,binsize,1:size(neuron.C,1),thresh,temp,[],[],countTimeThresh,small_velo);
