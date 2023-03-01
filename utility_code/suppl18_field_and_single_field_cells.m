function [firingrateAll,binTime,all_field_nums,single_field_idx]=suppl18_field_and_single_field_cells(foldername_all,behavName_all)

for j=1:length(foldername_all)
    load([foldername_all{j},'\','neuronIndividuals_new.mat'])
    for k=1:2
        load([foldername_all{j},'\',behavName_all{k}]);

        neuron=neuronIndividuals_new{k}.copy;
        behavpos=behav.position;
        behavtime=behav.time;
        maxbehavROI=behav.ROI;
        binsize=10;
        thresh=3*std(C_to_peakS(neuron.C),[],2);
        countTimeThresh=[0 inf];
        small_velo=10;

        [firingrateAll{j,k},~,~,binTime{j,k}] = calculatingCellSpatialForSingleData_040321(neuron,behavpos,behavtime,maxbehavROI,binsize,1:size(neuron.C,1),thresh,'S',[],[],countTimeThresh,small_velo);

        numFields=field_number_calculation(firingrateAll{j,k});

        all_field_nums{j,k}=numFields';
        single_field_idx{j,k}=find(numFields==1);
    end
end