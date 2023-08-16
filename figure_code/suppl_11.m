% suppl 11, cross correlation
Fs=15;
%% 1. calculate pearson correlation at border

% multiGeo
foldername_multiGeo={
    'D:\Xu_clusterting_paper_prep11_2020\arranged_final_data\Fig_1_3_tri_cir_sqr\M3411'	
    'D:\Xu_clusterting_paper_prep11_2020\arranged_final_data\Fig_1_3_tri_cir_sqr\M3412'	
    'D:\Xu_clusterting_paper_prep11_2020\arranged_final_data\Fig_1_3_tri_cir_sqr\M3421F'	
    'D:\Xu_clusterting_paper_prep11_2020\arranged_final_data\Fig_1_3_tri_cir_sqr\M3422F'	
    'D:\Xu_clusterting_paper_prep11_2020\arranged_final_data\Fig_1_3_tri_cir_sqr\M3424F'	
    'D:\Xu_clusterting_paper_prep11_2020\arranged_final_data\Fig_1_3_tri_cir_sqr\M3425F'	
    }

load('D:\Xu_clusterting_paper_prep11_2020\final_code\final_cluster_data\cluster_optimal_num\Fig1_3_multiGeo_clust_original.mat')

[all_corr_intra_multiGeo,all_corr_intra_border_multiGeo,all_corr_inter_border_multiGeo]=cluster_crossCorr2_border(group_ori_multiGeo,foldername_multiGeo)            

% Fig2
foldername_fig2={
    'D:\Xu_clusterting_paper_prep11_2020\final data\final_neuron_behav_data\Fig_2_cir_rec\M3411'		
    'D:\Xu_clusterting_paper_prep11_2020\final data\final_neuron_behav_data\Fig_2_cir_rec\M3412'		
    'D:\Xu_clusterting_paper_prep11_2020\final data\final_neuron_behav_data\Fig_2_cir_rec\M3413'		
    'D:\Xu_clusterting_paper_prep11_2020\final data\final_neuron_behav_data\Fig_2_cir_rec\M3414'		
    'D:\Xu_clusterting_paper_prep11_2020\final data\final_neuron_behav_data\Fig_2_cir_rec\M3415'		
    'D:\Xu_clusterting_paper_prep11_2020\final data\final_neuron_behav_data\Fig_2_cir_rec\M3421F'		
    'D:\Xu_clusterting_paper_prep11_2020\final data\final_neuron_behav_data\Fig_2_cir_rec\M3422F'		
    'D:\Xu_clusterting_paper_prep11_2020\final data\final_neuron_behav_data\Fig_2_cir_rec\M3423F'		
    'D:\Xu_clusterting_paper_prep11_2020\final data\final_neuron_behav_data\Fig_2_cir_rec\M3424F'		
    'D:\Xu_clusterting_paper_prep11_2020\final data\final_neuron_behav_data\Fig_2_cir_rec\M3425F'		
    'D:\Xu_clusterting_paper_prep11_2020\final data\final_neuron_behav_data\Fig_2_cir_rec\M3426F'		
    'D:\Xu_clusterting_paper_prep11_2020\final data\final_neuron_behav_data\Fig_2_cir_rec\M3427F'		
    }
load('D:\Xu_clusterting_paper_prep11_2020\final data\final_cluster_data\Fig2_cir_rec_clust_original.mat')

[all_corr_intra_Fig2,all_corr_intra_border_Fig2,all_corr_inter_border_Fig2]=cluster_crossCorr2_border(group_ori_fig2,foldername_fig2)            

% AI163
foldername_AI163={
    'D:\Xu_clusterting_paper_prep11_2020\final data\final_neuron_behav_data\Fig_3_barrier\101F\'
    'D:\Xu_clusterting_paper_prep11_2020\final data\final_neuron_behav_data\Fig_3_barrier\102F\'
    'D:\Xu_clusterting_paper_prep11_2020\final data\final_neuron_behav_data\Fig_3_barrier\103F\'
    'D:\Xu_clusterting_paper_prep11_2020\final data\final_neuron_behav_data\Fig_3_barrier\2833M\'
    'D:\Xu_clusterting_paper_prep11_2020\final data\final_neuron_behav_data\Fig_3_barrier\840F\'
}
load('D:\Xu_clusterting_paper_prep11_2020\final data\final_cluster_data\Fig3_barrier_clust_original.mat')

[all_corr_intra_AI163,all_corr_intra_border_AI163,all_corr_inter_border_AI163]=cluster_crossCorr2_border(group_ori_AI163,foldername_AI163)            

%% 3. illustration
% multiGeo
load('D:\Xu_clusterting_paper_prep11_2020\final_code\final_cluster_data\cluster_optimal_num\Fig1_3_multiGeo_clust_original.mat')
cluster_crossCorr_illustration2(group_ori_multiGeo,foldername_multiGeo,1,1,1,1,2)
cluster_crossCorr_trace_illustration(group_ori_multiGeo,foldername_multiGeo,1,1)

% Fig2
load('D:\Xu_clusterting_paper_prep11_2020\final data\final_cluster_data\Fig2_cir_rec_clust_original.mat')
cluster_crossCorr_illustration2(group_ori_fig2,foldername_fig2,8,1,11,3,2)
cluster_crossCorr_trace_illustration(group_ori_fig2,foldername_fig2,8,11)

% AI163
load('D:\Xu_clusterting_paper_prep11_2020\final data\final_cluster_data\Fig3_barrier_clust_original.mat')
cluster_crossCorr_illustration2(group_ori_AI163,foldername_AI163,2,1,2,1,2)
cluster_crossCorr_trace_illustration(group_ori_AI163,foldername_AI163,2,2)

%% 4. pearson corr illustration 
colorClusters_all=colormap(lines);
% multiGeo

subplot(131)

h1=cdfplot(cell2mat(all_corr_intra_multiGeo(:,3)));
hold on;
h2=cdfplot(cell2mat(all_corr_intra_border_multiGeo(:,3)));
h3=cdfplot(cell2mat(all_corr_inter_border_multiGeo(:,3)));
set(h1,'color',colorClusters_all(4,:));
set(h2,'color',colorClusters_all(5,:));
set(h3,'color',colorClusters_all(6,:));

p1=infer_cdf_loc(cell2mat(all_corr_intra_multiGeo(:,3)),mean(cell2mat(all_corr_intra_multiGeo(:,3))));
p2=infer_cdf_loc(cell2mat(all_corr_intra_border_multiGeo(:,3)),mean(cell2mat(all_corr_intra_border_multiGeo(:,3))));
p3=infer_cdf_loc(cell2mat(all_corr_inter_border_multiGeo(:,3)),mean(cell2mat(all_corr_inter_border_multiGeo(:,3))));

plot(p1(1),p1(2),'.','color',colorClusters_all(4,:),'MarkerSize',20)
plot(p2(1),p2(2),'.','color',colorClusters_all(5,:),'MarkerSize',20)
plot(p3(1),p3(2),'.','color',colorClusters_all(6,:),'MarkerSize',20)

[v,p]=intra_inter_statistic({all_corr_intra_multiGeo(:,3),all_corr_intra_border_multiGeo(:,3),all_corr_inter_border_multiGeo(:,3)},{'intra-cluster','border intra-cluster','border inter-cluster'});

% multiGeo, per mice ver
[avg_corr_mice_intra_multiGeo]=corr_border_processing_perMice(all_corr_intra_multiGeo(:,3));
[avg_corr_mice_intra_border_multiGeo]=corr_border_processing_perMice(all_corr_intra_border_multiGeo(:,3));
[avg_corr_mice_inter_border_multiGeo]=corr_border_processing_perMice(all_corr_inter_border_multiGeo(:,3));

% fig2

subplot(132)

h1=cdfplot(cell2mat(all_corr_intra_Fig2(:,1)));
hold on;
h2=cdfplot(cell2mat(all_corr_intra_border_Fig2(:,1)));
h3=cdfplot(cell2mat(all_corr_inter_border_Fig2(:,1)));
set(h1,'color',colorClusters_all(4,:));
set(h2,'color',colorClusters_all(5,:));
set(h3,'color',colorClusters_all(6,:));

p1=infer_cdf_loc(cell2mat(all_corr_intra_Fig2(:,1)),mean(cell2mat(all_corr_intra_Fig2(:,1))));
p2=infer_cdf_loc(cell2mat(all_corr_intra_border_Fig2(:,1)),mean(cell2mat(all_corr_intra_border_Fig2(:,1))));
p3=infer_cdf_loc(cell2mat(all_corr_inter_border_Fig2(:,1)),mean(cell2mat(all_corr_inter_border_Fig2(:,1))));

plot(p1(1),p1(2),'.','color',colorClusters_all(4,:),'MarkerSize',20)
plot(p2(1),p2(2),'.','color',colorClusters_all(5,:),'MarkerSize',20)
plot(p3(1),p3(2),'.','color',colorClusters_all(6,:),'MarkerSize',20)

[v,p]=intra_inter_statistic({all_corr_intra_Fig2(:,1),all_corr_intra_border_Fig2(:,1),all_corr_inter_border_Fig2(:,1)},{'intra-cluster','border intra-cluster','border inter-cluster'});

% AI163

subplot(133)

h1=cdfplot(cell2mat(all_corr_intra_AI163(:,1)));
hold on;
h2=cdfplot(cell2mat(all_corr_intra_border_AI163(:,1)));
h3=cdfplot(cell2mat(all_corr_inter_border_AI163(:,1)));
set(h1,'color',colorClusters_all(4,:));
set(h2,'color',colorClusters_all(5,:));
set(h3,'color',colorClusters_all(6,:));

p1=infer_cdf_loc(cell2mat(all_corr_intra_AI163(:,1)),mean(cell2mat(all_corr_intra_AI163(:,1))));
p2=infer_cdf_loc(cell2mat(all_corr_intra_border_AI163(:,1)),mean(cell2mat(all_corr_intra_border_AI163(:,1))));
p3=infer_cdf_loc(cell2mat(all_corr_inter_border_AI163(:,1)),mean(cell2mat(all_corr_inter_border_AI163(:,1))));

plot(p1(1),p1(2),'.','color',colorClusters_all(4,:),'MarkerSize',20)
plot(p2(1),p2(2),'.','color',colorClusters_all(5,:),'MarkerSize',20)
plot(p3(1),p3(2),'.','color',colorClusters_all(6,:),'MarkerSize',20)

[v,p]=intra_inter_statistic({all_corr_intra_AI163(:,1),all_corr_intra_border_AI163(:,1),all_corr_inter_border_AI163(:,1)},{'intra-cluster','border intra-cluster','border inter-cluster'});

