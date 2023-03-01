function neuron_cat=neuronIndividuals_new_concatenate(neuronIndividuals_new)

    for j1=1:length(neuronIndividuals_new)
        neuron_cat=Sources2D;
        neuron_cat.A=neuronIndividuals_new{j1}.A;
        neuron_cat.Cn=neuronIndividuals_new{j1}.Cn;
        neuron_cat.Coor=neuronIndividuals_new{j1}.Coor;
        neuron_cat.centroid=neuronIndividuals_new{j1}.centroid;
        
        neuron_cat.C=[neuron_cat.C;neuronIndividuals_new{j1}.C];
        neuron_cat.S=[neuron_cat.S;neuronIndividuals_new{j1}.S];
        neuron_cat.C_raw=[neuron_cat.S;neuronIndividuals_new{j1}.C_raw];
        
        timeDiff=round(mean(diff(neuron_cat.time)));
        if isnan(timeDiff)
            timeDiff=0;
        end
        neuron_cat.time=[neuron_cat.time;neuronIndividuals_new{j1}.time+timeDiff];
    end