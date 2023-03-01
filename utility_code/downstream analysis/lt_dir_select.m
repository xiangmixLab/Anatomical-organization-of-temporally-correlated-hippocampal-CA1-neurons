function [all_pc_t,all_ratemaps,all_cts]=lt_dir_select(all_neuron_simp_laps1,all_behav_laps1,all_pc_dir,all_neuron_num,binsize,infoscore_type,dir_sign)

if dir_sign==1  
    all_pc_dirt=all_pc_dir{1};
else
    all_pc_dirt=all_pc_dir{2};
end
    
for i=1:size(all_neuron_simp_laps1,1)
    for j=1:size(all_neuron_simp_laps1,2)
        [fr_all_dir{i,j},ct_all_dir{i,j},ctime_all_dir{i,j}] = calculatinglinearTrackRateMaps_laps(all_neuron_simp_laps1{i,j},all_behav_laps1{i,j},binsize,'S',all_neuron_num(i));
        
        all_ratemaps{i,j}=fr_all_dir{i,j}{dir_sign}{1};
        all_cts{i,j}=ctime_all_dir{i,j}{dir_sign}{1};
        
        if size(all_cts{i,j},2)>size(all_cts{i,j},1) % hor track
            for k=1:length(all_ratemaps{i,j})
                all_ratemaps{i,j}{k}=nanmean(all_ratemaps{i,j}{k},1);
            end
            all_cts{i,j}=nanmean(all_cts{i,j},1);
        else % vet track
           for k=1:length(all_ratemaps{i,j})
                all_ratemaps{i,j}{k}=nanmean(all_ratemaps{i,j}{k},2);
            end
            all_cts{i,j}=nanmean(all_cts{i,j},2);    
        end
    end
end

all_pc_t={};
for i=1:size(all_pc_dirt,1)
    for j=1:size(all_pc_dirt,2)
        all_pc_t{i,j}=all_pc_dirt{i,j}{infoscore_type};
    end
end