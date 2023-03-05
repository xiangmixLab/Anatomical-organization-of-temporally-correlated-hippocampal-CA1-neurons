% cluster spatial connected region analysis
% use DBSCAN to isolate locally compacted neuron populations
% then determine the boundary cells and generate the region
function [boundary_all,region_all,region_nnum_all,nd,boundary_idx_all,A_color_test]=cluster_region_cal_022422(group1,neuronIndividuals_new,neighboor_dis_given)

d1=neuronIndividuals_new{1}.imageSize(1);
d2=neuronIndividuals_new{1}.imageSize(2);

centroidd=neuron_centroid_calculation(neuronIndividuals_new{1},[d1,d2]);

gp1=group1(group1~=-1);

%% 1: determine neighborhood distance, across all neurons
neighboor_dis_all=squareform(pdist(centroidd));
neighboor_dis_all(neighboor_dis_all==0)=inf;
neighboor_dis=quantile(min(neighboor_dis_all,[],2),0.95);

if ~isempty(neighboor_dis_given)&&neighboor_dis<=neighboor_dis_given
    neighboor_dis=neighboor_dis_given;
end

%% 2. use DBSCAN to calculate spatially clustered neuron subgroups
boundary_all={};
boundary_idx_all={};

region_nnum_all={};

for gp=1:length(unique(gp1))

    cen=round(centroidd(group1==gp,:));
    nd=neighboor_dis;

    [IDX, isnoise]=DBSCAN_modified(cen,nd,4,'dis');
    
    uIDX=unique(IDX);
    uIDX(uIDX==0)=[];
    
    for i=1:length(unique(IDX))
        cen_t=[cen(IDX==i,1),cen(IDX==i,2)];
        k=boundary(cen_t(:,1),cen_t(:,2),0.9);
        boundary_all{gp}{i}=[cen_t(k,1),cen_t(k,2)];
        region_nnum_all{gp}{i}=sum(IDX==i);
        
        curClustIdx=find(IDX==i);
        globalIdx=find(group1==gp);
        boundary_idx_all{gp}{i}=globalIdx(curClustIdx(k));
    end

end

%% 3: fill boundary to get regions
region_all={};
for gp=1:length(unique(gp1))
    for i=1:length(boundary_all{gp})
        boundd=boundary_all{gp}{i};
        Ai=poly2mask(boundd(:,1),boundd(:,2),d1,d2); % or the regions are flipped when showing
        region_all{gp}{i}=Ai;
    end
end


% test
colorClusters_all=distinguishable_colors(10);
A_color_test = cluster_spatial_footprint_colormap(neuronIndividuals_new,d1,d2,colorClusters_all,group1,0.8); 

% region detection test
% figure;axis ij;hold on;
% color_cluster=distinguishable_colors(10);
% for gp=1:length(unique(gp1))
%     plot(centroidd(group1==gp,1),centroidd(group1==gp,2),'.','markerSize',5,'color',color_cluster(gp,:));
% end
% 
% color_cluster=distinguishable_colors(10);
% for gp=1:length(unique(gp1))
%     for i=1:length(boundary_all{gp})
%         plot(boundary_all{gp}{i}(:,1),boundary_all{gp}{i}(:,2),'color',color_cluster(gp,:));
%     end
% end