function [neuron_combined,behav_combined]=OF_trial_combine_data(neuron,behav,del_idx)
    
    neuron_combined.C=[];
    neuron_combined.C_raw=[];
    neuron_combined.S=[];
    neuron_combined.time=[];
    
    behav_combined.position=[];
    behav_combined.time=[];
    behav_combined.ROI=[];
    
    if ~isempty(neuron)
        timet=[0];
        for j=1:length(neuron)
           neuron_combined.C=[neuron_combined.C,neuron{j}.C(setdiff(1:size(neuron{j}.C,1),del_idx),:)];
           neuron_combined.C_raw=[neuron_combined.C_raw,neuron{j}.C_raw(setdiff(1:size(neuron{j}.C,1),del_idx),:)];
           timet=[timet;neuron{j}.time+timet(end)+mean(diff(neuron{j}.time))];
        end

        neuron_combined.S=C_to_peakS(neuron_combined.C);
        neuron_combined.time=timet(2:end);
    end
    
    if ~isempty(behav)
        timet=[0];
        ROI_all=[];
        for j=1:length(behav)
           behav_combined.position=[behav_combined.position;behav{j}.position];
           timet=[timet;behav{j}.time+timet(end)+mean(diff(behav{j}.time))];
           ROI_all(j,:)=behav{j}.ROI;
        end
        
        behav_combined.time=timet(2:end);
        behav_combined.ROI=[0 0 mean(ROI_all(:,3:4),1)];
    end
    