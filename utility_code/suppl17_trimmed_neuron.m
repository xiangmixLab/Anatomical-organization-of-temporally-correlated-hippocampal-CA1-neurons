function suppl17_trimmed_neuron(dpath,preserve_idxs)

for i=1:length(dpath)
    load([dpath{i},'\','neuronIndividuals_new.mat'])
    load([dpath{i},'\','neuronIndividuals_new_tuncat.mat'])
    
    for j=1:length(neuronIndividuals_new)
        neuronIndividuals_new{j}.delete(find(preserve_idxs{i,j}==0));
        neuronIndividuals_new_tuncat{j}.delete(find(preserve_idxs{i,j}==0));
    end
    
    save([dpath{i},'\','neuronIndividuals_new_cleaned.mat'],'neuronIndividuals_new','-v7.3');
    save([dpath{i},'\','neuronIndividuals_new_tuncat_cleaned.mat'],'neuronIndividuals_new_tuncat','-v7.3');
end