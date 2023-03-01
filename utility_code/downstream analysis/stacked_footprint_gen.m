%% stack footprint is a 3d stack in which each frame contains the footprint of 1 neuron
function all_footprints=stacked_footprint_gen(neuron)

all_footprints=zeros(neuron.imageSize(1),neuron.imageSize(2),size(neuron.A,2));
for i=1:size(neuron.A,2)
    curr=reshape(neuron.A(:,i),neuron.imageSize(1),neuron.imageSize(2));
    curr=curr>max(curr(:))*0.1;
    all_footprints(:,:,i)=curr;
end

all_footprints=logical(all_footprints);
    