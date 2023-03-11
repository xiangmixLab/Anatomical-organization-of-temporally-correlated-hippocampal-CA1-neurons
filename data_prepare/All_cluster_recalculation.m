%% calculate clusters for all experiments used in study
run('neuron_data_info.m');
%% 1. cluster calculation: each trial has its own cluster number
% multiGeo
[group_multiGeo,group_ori_multiGeo]=experiment_cluster_cal(foldername_multiGeo,6,6,[]);

% cir-rec
[group_fig2,group_ori_fig2]=experiment_cluster_cal(foldername_fig2,12,2,[]);

% barrier
[group_AI163,group_ori_AI163]=experiment_cluster_cal(foldername_AI163,5,9,[]);

% linear track
[group_LT,group_ori_LT]=experiment_cluster_cal(foldername_LT,6,3,[]);

% immobile
[group_imb,group_ori_imb]=experiment_cluster_cal(foldername_immobile,3,2,[]);

%% 2. cluster calculation: cluster number 2-10

% multiGeo
[group_multiGeo_k,group_ori_multiGeo_k]=experiment_cluster_cal_2_to_10(foldername_multiGeo,6,6);

% cir-rec
[group_fig2_k,group_ori_fig2_k]=experiment_cluster_cal_2_to_10(foldername_fig2,12,2);

% barrier
[group_AI163_k,group_ori_AI163_k]=experiment_cluster_cal_2_to_10(foldername_AI163,5,9);

% linear track
[group_LT_k,group_ori_LT_k]=experiment_cluster_cal_2_to_10(foldername_LT,6,3);

% immobile
[group_imb_k,group_ori_imb_k]=experiment_cluster_cal_2_to_10(foldername_immobile,3,2);

%% 3. cluster calculation for figure 3's overlap: (multiGeometry have 5 clusters, linear track have 4 clusters, barrier have 4 clusters, as indicated in method)

% multiGeo
[group_multiGeo_sameNum_allMice,group_ori_multiGeo_sameNum_allMice]=experiment_cluster_cal(foldername_multiGeo,6,6,5);

% linear track
[group_LT_sameNum_allMice,group_ori_LT_sameNum_allMice]=experiment_cluster_cal(foldername_LT,6,3,4);

% barrier
[group_AI163_sameNum_allMice,group_ori_AI163_sameNum_allMice]=experiment_cluster_cal(foldername_AI163,5,9,4);


%% 4. cluster calculation for figure 3's overlap, first and second half

% multiGeo
[group_multiGeo_half,group_ori_multiGeo_half]=experiment_cluster_cal_half(foldername_multiGeo,6,6,5);

% linear track
[group_LT_half,group_ori_LT_half]=experiment_cluster_cal_half(foldername_LT,6,3,4);

% barrier
[group_AI163_half,group_ori_AI163_half]=experiment_cluster_cal_half(foldername_AI163,5,9,4);

%% 5. cluster size calculation, cluster num 2-10

% multiGeo
[A_color_multiGeo_k,A_color_region_multiGeo_k,avg_region_multiGeo_k,avg_region_shuf_multiGeo_k]=experiment_cluster_region_cal_2_to_10(group_ori_multiGeo_k,foldername_multiGeo);

% cir-rec
[A_color_fig2_k,A_color_region_fig2_k,avg_region_fig2_k,avg_region_shuf_fig2_k]=experiment_cluster_region_cal_2_to_10(group_ori_fig2_k,foldername_fig2);

% barrier
[A_color_AI163_k,A_color_region_AI163_k,avg_region_AI163_k,avg_region_shuf_AI163_k]=experiment_cluster_region_cal_2_to_10(group_ori_AI163_k,foldername_AI163);

% linear track
[A_color_LT_k,A_color_region_LT_k,avg_region_LT_k,avg_region_shuf_LT_k]=experiment_cluster_region_cal_2_to_10(group_ori_LT_k,foldername_LT);

% immobile
[A_color_imb_k,A_color_region_imb_k,avg_region_imb_k,avg_region_shuf_imb_k]=experiment_cluster_region_cal_2_to_10(group_ori_imb_k,foldername_immobile);

