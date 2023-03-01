function [esm_ratemap_all,esm_ratemap_clust_2to1,esm_ratemap_clust_1to2,esm_ratemap_nonsmooth_all,A_color_2to1,A_color_1to2]=determineReorganizedEnsembleRatemap_072922(group1,group2,group12,group21,ratemap1,ratemap2)

ugp1=unique(group1);
ugp2=unique(group2);
ugp1(ugp1==-1)=[];
ugp2(ugp2==-1)=[];

esm_ratemap_all=cell(2,max(max(ugp1),max(ugp2)));
esm_ratemap_nonsmooth_all=cell(2,max(max(ugp1),max(ugp2)));

[esm_ratemap_all(1,ugp1),esm_ratemap_nonsmooth_all(1,ugp1)]=clust_ensemble_ratemap(group1,ratemap1);
[esm_ratemap_all(2,ugp2),esm_ratemap_nonsmooth_all(2,ugp2)]=clust_ensemble_ratemap(group2,ratemap2);

esm_ratemap_clust_2to1=cell(max(ugp2),max(ugp1)); % check how cluster 1 is represented by cluster 2 (how cluster 1's ensemble ratemap fractured/expanded)
A_color_2to1={};
for i=1:size(group21,1)
    idxs=group21{i,3};
    
    idxs_mat=[];
    gp_idxs_mat=[];
    for j=1:length(idxs)
        idxs_mat=[idxs_mat;idxs{j}];
        gp_idxs_mat=[gp_idxs_mat;ones(length(idxs{j}),1)*j];
    end
    esm_ratemap_clust_2to1(i,unique(gp_idxs_mat))=clust_ensemble_ratemap(gp_idxs_mat,ratemap2(idxs_mat));
%     A_color_2to1{i}=cluster_spatial_footprint_colormap(n1,n1.imageSize(1),n1.imageSize(2),colorClusters_all,group1,0.75);

end

esm_ratemap_clust_1to2=cell(max(ugp1),max(ugp2)); % check how cluster 1 is represented by cluster 2 (how cluster 1's ensemble ratemap fractured/expanded)
A_color_1to2={};
for i=1:size(group12,1)
    idxs=group12{i,3};
    
    idxs_mat=[];
    gp_idxs_mat=[];
    for j=1:length(idxs)
        idxs_mat=[idxs_mat;idxs{j}];
        gp_idxs_mat=[gp_idxs_mat;ones(length(idxs{j}),1)*j];
    end
    esm_ratemap_clust_1to2(i,unique(gp_idxs_mat))=clust_ensemble_ratemap(gp_idxs_mat,ratemap1(idxs_mat));
end
