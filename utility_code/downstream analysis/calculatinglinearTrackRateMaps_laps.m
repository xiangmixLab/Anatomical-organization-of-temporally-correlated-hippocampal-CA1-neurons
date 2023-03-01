function [fr_all_dir,ct_all_dir,ctime_all_dir,fr_all,ctime_all] = calculatinglinearTrackRateMaps_laps(ndat,bdat,binsize,temp,all_neuron_num)


for k1=1:size(ndat,1)
    for k2=1:size(ndat,2)
        if ~isempty(ndat{k1,k2})
            [fr_all{1,1}{k1,k2},ct_all{1,1}{k1,k2},~,ctime_all{1,1}{k1,k2}]=calculatingCellSpatialForSingleData_040321(ndat{k1,k2},bdat{k1,k2}.position,bdat{k1,k2}.time,bdat{k1,k2}.ROI,binsize,1:size(ndat{k1,k2}.C,1),3*std(ndat{k1,k2}.S,[],2),temp,[],[],[0 inf],0);
        end
    end
end


[fr_all_dir,ct_all_dir,ctime_all_dir]=laps_avg_fr_ct(fr_all,ct_all,ctime_all,all_neuron_num);
