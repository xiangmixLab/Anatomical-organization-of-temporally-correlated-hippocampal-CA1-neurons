%DBSCAN region properities for different num of clusters
function [A_color,A_color_region,avg_region,max_region,avg_region_nnum,max_region_nnum,avg_region_shuffled,max_region_shuffled,avg_region_nnum_shuffled,max_region_nnum_shuffled,nd,A_color_shuf,A_color_region_shuf]=DBSCAN_region_quantify_noShuf_041522(group1,neuronIndividuals_new,neighboor_dis_given)

colorClusters_all=distinguishable_colors(10);
shuffles=100;
if isempty(neighboor_dis_given)
%     neighboor_dis_given=17.5; % 35um, Dombeck et al. 2010
    neighboor_dis_given=-1; % use 95% of neighborhood dis
end

%% 1. calculate spatial compact regions with DBSCAN
[boundaryy,regionn,region_nnum,nd]=cluster_region_cal_022422(group1,neuronIndividuals_new,neighboor_dis_given);

% shuffled case
group_rand={};
for k=1:shuffles
    group_rand{k}= group1(randperm(length(group1)));
%     A_color_rand{k} = cluster_spatial_footprint_colormap(neuronIndividuals_new,240,376,colorClusters_all,group_rand,0.8); 
    [boundaryy_shuffled{k},regionn_shuffled{k},region_nnum_shuffled{k}]=cluster_region_cal_022422(group_rand{k},neuronIndividuals_new,neighboor_dis_given);      
end

%% 2. calculate region statistics
imgSize=neuronIndividuals_new{1}.imageSize;

avg_region_shuffled=[];
max_region_shuffled=[];
area_thresh_all_shuffled={};

[avg_region,max_region,avg_region_nnum,max_region_nnum,max_regions,area_thresh_all]=DBSCAN_region_valueCal_022422(regionn,region_nnum,imgSize);
for k=1:shuffles
    [avg_region_shuffled(k),max_region_shuffled(k),avg_region_nnum_shuffled(k),max_region_nnum_shuffled(k),area_thresh_all_shuffled{k}]=DBSCAN_region_valueCal_022422(regionn_shuffled{k},region_nnum_shuffled{k},imgSize);
end

%% 3. footprint and regions, normal clustering
[A_color,A_color_region]=DBSCAN_footprint_regions(neuronIndividuals_new,imgSize,group1,regionn,area_thresh_all);

%% 4. footprint and regions, shuffle clustering
A_color_shuf={};
A_color_region_shuf={};
% parfor k=1:shuffles
%     [A_color_shuf{k},A_color_region_shuf{k}]=DBSCAN_footprint_regions(neuronIndividuals_new,imgSize,group_rand{k},regionn_shuffled{k},area_thresh_all);
% end
