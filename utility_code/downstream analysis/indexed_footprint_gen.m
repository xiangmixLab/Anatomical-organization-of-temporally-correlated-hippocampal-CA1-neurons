function all_footprints=indexed_footprint_gen(neuron)

all_footprints=zeros(neuron.imageSize(1),neuron.imageSize(2));
for i=1:size(neuron.A,2)
    curr=reshape(neuron.A(:,i),neuron.imageSize(1),neuron.imageSize(2));
    curr=curr>max(curr(:))*0.1;
    all_footprints=all_footprints+curr*i;
end
    