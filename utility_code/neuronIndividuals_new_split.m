function neuronIndividuals=neuronIndividuals_new_split(neuron,frames)

startt=1;
startt1=1;
neuronIndividuals={};
for i=1:2*floor(size(neuron.C,2)/frames)-1
    neuronIndividuals{i}=neuron;
    neuronIndividuals{i}.C=neuronIndividuals{i}.C(:,startt:startt+frames-1);
    neuronIndividuals{i}.S=neuronIndividuals{i}.S(:,startt:startt+frames-1);
    neuronIndividuals{i}.C_df=[];
    neuronIndividuals{i}.trace=[];
    neuronIndividuals{i}.C_raw=neuronIndividuals{i}.C_raw(:,startt:startt+frames-1);
    neuronIndividuals{i}.time=neuronIndividuals{i}.time(startt1:startt1+2*frames-1);

    startt=startt+frames/2;
    startt1=startt1+2*frames/2;
end
