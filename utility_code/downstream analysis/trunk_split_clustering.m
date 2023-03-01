function [group_trunk,Acolor_trunk,Acolor_region_trunk,region_size_trunk,avg_region_size]=trunk_split_clustering(neuron0,trunk_size,gp_num,imageSize)

neuron_split={};
timeleng=size(neuron0.C,2);
neuron0.imageSize=imageSize;
[~,centroid]=spatial_footprint_calculation(neuron0);
for i=1:timeleng-trunk_size
    neuron_split{i}.A=full(neuron0.A);
    neuron_split{i}.C=neuron0.C(:,i:i+trunk_size-1);
    neuron_split{i}.S=neuron0.C(:,i:i+trunk_size-1);
    neuron_split{i}.imageSize=imageSize;
    
    neuron_split{i}.centroid=centroid;
end



group_trunk={};
Acolor_trunk={};
Acolor_region_trunk={};
region_size_trunk=[];

parfor i=1:length(neuron_split)
    group_trunk{i}=cluster_determine_by_suoqin_NMF_firstPeakCoph(neuron_split{i},100,10,gp_num);
    [Acolor_trunk{i},Acolor_region_trunk{i},region_size_trunk(i)]=DBSCAN_region_quantify_func_single_group_no_plot(group_trunk{i},neuron_split(i));
end

avg_region_size=nanmean(region_size_trunk);    

