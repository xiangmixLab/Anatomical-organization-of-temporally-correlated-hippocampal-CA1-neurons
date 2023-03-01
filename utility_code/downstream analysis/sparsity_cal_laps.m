function spr=sparsity_cal_laps(all_neuron,all_behav,binsize,temp,all_neuron_num,dir_sign)

[fr_all_dir,ct_all_dir,ctime_all_dir] = calculatinglinearTrackRateMaps_laps(all_neuron,all_behav,temp,binsize,all_neuron_num);

firingrateAllt=fr_all_dir{dir_sign}{1}; % dir_sign: 1: left 2: right
countTime=ctime_all_dir{dir_sign}{1}; % 

Po=countTime/sum(countTime(:));
spr=[];
for j=1:length(firingrateAllt)

    if ~isempty(firingrateAllt{j}) && sum(firingrateAllt{j}(:))~=0
%         firingrateAllt{j}=nansum(firingrateAllt{j},1);
%         Po=nansum(countTime,1)/sum(countTime(:));
        spr(j)=sum(sum(firingrateAllt{j}.*Po))^2/sum(sum(firingrateAllt{j}.^2.*Po));
    else
        spr(j)=nan;
    end
end
