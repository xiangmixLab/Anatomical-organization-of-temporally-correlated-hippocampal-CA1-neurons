%% HDAC AD processing script

function neuronIndividuals_new_generation_030521(foldernamet,timestampname,mscamid,behavcamid,mrange)

    %% start
    for ikk=mrange  

        cd(foldernamet{ikk})
        
        % read data
        [neuron,neuronfilename]=read_neuron_data(ikk,foldernamet,1);

        % change neuron.S
        neuron.S=C_to_peakS(neuron.C);

        % split neuron data
        neuronIndividuals_new = splittingMUltiConditionNeuronData_adapted_automatic_032121(neuron,neuronfilename,mscamid{ikk},behavcamid{ikk},timestampname); 
        
        % save
        save('neuronIndividuals_new.mat','neuronIndividuals_new','-v7.3');
    end