function asc=spatial_coherence_cal_laps(all_neuron,all_behav,binsize,temp,all_neuron_num,dir_sign)

[fr_all_dir,ct_all_dir,ctime_all_dir] = calculatinglinearTrackRateMaps_laps(all_neuron,all_behav,temp,binsize,all_neuron_num);

firingrateAllt=fr_all_dir{dir_sign}{1}; % dir_sign: 1: left 2: right

asc=[];
for k=1:length(firingrateAllt)
    if ~isempty(firingrateAllt{k})
%         firingrateAllt{k}=nansum(firingrateAllt{k},1);
        [~,asc(k,1)]=spatial_coherence(firingrateAllt{k});
    else
        asc(k,1)=nan;
    end
end