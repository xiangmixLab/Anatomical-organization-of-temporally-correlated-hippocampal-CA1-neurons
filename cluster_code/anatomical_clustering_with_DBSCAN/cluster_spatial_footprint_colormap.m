function A_color=cluster_spatial_footprint_colormap(neuronIndividuals_new,d1,d2,colorClusters_all,group,thresh)
A_color=ones(d1*d2,3);
for celll=1:size(neuronIndividuals_new{1}.A,2)
    Ai=reshape(neuronIndividuals_new{1}.A(:,celll),d1,d2);
    Ai(Ai<thresh*max(Ai(:)))=0;
    Ai=logical(Ai);
    Ai=full(Ai);
    Ai=bwareaopen(Ai,9);
    se=strel('disk',2);
    Ai=imclose(Ai,se);
    ind=find(Ai>0);%the pixels at which a neuron shows up
    if group(celll)~=-1
        A_color(ind,:)=repmat(colorClusters_all(group(celll),:),length(ind),1);
    else
        A_color(ind,:)=repmat([0.9,0.9,0.9],length(ind),1);
    end
end
A_color=reshape(A_color,d1,d2,3);