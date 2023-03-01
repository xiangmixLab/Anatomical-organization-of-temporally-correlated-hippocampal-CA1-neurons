function [fr_all_laps,ct_all_laps,ctime_all_laps]=calculate_laps_ratemaps(all_neuron_simp_laps,all_behav_laps,all_behav,binsize,smallvelo)
fr_all_laps={};
ct_all_laps={};
ctime_all_laps={};
for tk=1:size(all_neuron_simp_laps,1)
    for i=1:size(all_neuron_simp_laps,2)
        
        for j1=1:size(all_neuron_simp_laps{tk,i}.C,1)
            for j2=1:size(all_neuron_simp_laps{tk,i}.C,2)
                if ~isempty(all_neuron_simp_laps{tk,i}.C{j1,j2})
                    ansl_temp.C=all_neuron_simp_laps{tk,i}.C{j1,j2};
                    ansl_temp.S=all_neuron_simp_laps{tk,i}.S{j1,j2};
                    ansl_temp.time=all_neuron_simp_laps{tk,i}.time{j1,j2};
                    [fr_all_laps{tk,i}{j1,j2},ct_all_laps{tk,i}{j1,j2},~,ctime_all_laps{tk,i}{j1,j2}]=calculatingCellSpatialForSingleData_Suoqin(ansl_temp,all_behav_laps{tk,i}.position{j1,j2},all_behav_laps{tk,i}.time{j1,j2},all_behav{tk,i}.ROI,binsize,1:size(ansl_temp.C,1),0.1*max(ansl_temp.C,[],2),'S',[],[],[0 inf],smallvelo);

                end
            end
        end
    end
end
