function [neuron1,neuron_del]=CNMFE_manual_intervene(neuron)
    
    neuron1=neuron.copy;
    show_merge = true;
    neuron1.orderROIs('snr');   % order neurons in different ways {'snr', 'decay_time', 'mean', 'circularity'}
    ind_del=neuron1.viewNeurons([], neuron1.C_raw);

    % merge closeby neurons
%     dmin_only = 2;
%     [~,newIDs]=neuron1.merge_close_neighbors(true, dmin_only);

    % delete neurons
    tags = neuron1.tag_neurons_parallel();  % find neurons with fewer nonzero pixels than min_pixel and silent calcium transients
    ids = find(tags>0); 
    if ~isempty(ids)
        neuron1.viewNeurons(ids, neuron1.C_raw);
    end
    
    neuron_del=neuron.copy;
    neuron_del.delete(setdiff(1:size(neuron.C,1),ind_del));