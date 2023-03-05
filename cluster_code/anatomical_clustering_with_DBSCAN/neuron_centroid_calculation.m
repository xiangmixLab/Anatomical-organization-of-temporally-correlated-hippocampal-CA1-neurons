function centroidd=neuron_centroid_calculation(neuron,dsize)

try
neuron.A=full(neuron.A);
catch
end

centroidd=[];
for i=1:size(neuron.A,2)
    At=neuron.A(:,i);
    At=reshape(At,dsize(1),dsize(2));
    At(At<0.5*max(At(:)))=0;
    At=logical(At);
    stats=regionprops(At);
    centroidd(i,:)=stats.Centroid;
end