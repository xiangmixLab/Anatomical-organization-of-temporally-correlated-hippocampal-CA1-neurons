function [A_color,A_color_region]=DBSCAN_footprint_regions(neuronIndividuals_new,imgSize,group1,regionn,area_thresh_all)

colorClusters_all=distinguishable_colors(10);
%% 3-1.inneed: footprints
A_color=cluster_spatial_footprint_colormap(neuronIndividuals_new,imgSize(1),imgSize(2),colorClusters_all,group1,0.75);

%% 3-2.inneed: regions
A_color_region=ones(imgSize(1)*imgSize(2),3)*255;
for i=1:length(regionn)
    for j=1:length(regionn{i})
        if ~isempty(regionn{i}{j})
            region_t=regionn{i}{j}>0;
            region_t=bwareaopen(region_t,round(area_thresh_all(i)));
            ind=find(region_t>0);%the pixels at which a neuron shows up
            A_color_region(ind,:)=repmat(colorClusters_all(i,:),length(ind),1);
        end
    end
end
A_color_region=reshape(A_color_region,imgSize(1),imgSize(2),3);
