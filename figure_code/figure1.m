run('D:\Xu_clusterting_paper_prep11_2020\final_code\data_prepare\neuron_data_info.m')
all_colors=distinguishable_colors(10);
load(['D:\Xu_clusterting_paper_prep11_2020\final_code\final_cluster_data\cluster_optimal_num\Fig1_3_multiGeo_clust_original.mat'])

%% A. behav
load([foldername_multiGeo{6},'\','behav.mat']);
plot(behavIndividuals{2}.position(:,1),behavIndividuals{2}.position(:,2),'color',[0.5 0.5 0.5],'lineWidth',2);

%% B: footprint
load([foldername_multiGeo{5},'\','neuronIndividuals_new.mat'])
footprint=spatial_footprint_calculation(neuronIndividuals_new{1});
imagesc(footprint)

%% C: traces, organized by clust
load([foldername_multiGeo{6},'\','neuronIndividuals_new.mat'])

group_curr=group_ori_multiGeo{6,2};
clustered_timeSeries2(neuronIndividuals_new{2}.C,group_curr)

%% D1: magify trace: performed in Adobe illustrator, direct magnify the circled segment in C
dataC3 = zeros(size(neuronIndividuals_new{2}.C,1),751);
segment={[750:1500],[1800:2550],[2000:2750],[6300:7050]};
for j = 1:length(unique(group_curr))
    dataC3(group_curr == j,:) = neuronIndividuals_new{2}.C(group_curr == j,segment{j});
end
clustered_timeSeries2(dataC3,group_curr)

    
%% D2: trajectories during magnified period
load([foldername_multiGeo{6},'\','behav.mat']);

figure;
behavIndividuals{2}.position=behavIndividuals{2}.position*280/299;
for i=1:length(segment)
    plot(behavIndividuals{2}.position(segment{i},1),behavIndividuals{2}.position(segment{i},2),'color',all_colors(i,:),'lineWidth',2);
    traj_leng(i)=trajectory_length(behavIndividuals{2}.position(segment{i},:));
    hold on;
end
cen=mean(behavIndividuals{2}.position,1)-[12 1]; % mean do not directly represent the real center, eyeball adjustment
viscircles(cen,140);

%% E: correlation map
perm=[];
for j = 1:length(unique(group_curr))
    perm = [perm;find(group_curr == j)];
end
CM22=squareform(1-pdist(neuronIndividuals_new{2}.C,'correlation'));
displayHeatmap_general(CM22,CM22,perm,all_colors,max(group_curr),group_curr,1);

%% F: Intra, inter cluster correlation, wiht shuffle, circle, use cumulative distribution
for i=1:6
    load([foldername_multiGeo{i},'\','neuronIndividuals_new.mat'])
    group_curr=group_ori_multiGeo{i,2};
    [intra_all{i},inter_all{i},~,~,intra_shuffle_all{i},~]=intra_inter_cluster_corr_dis(neuronIndividuals_new,group_curr,2,'corr');
end

cdfplot(cell2mat(intra_all')); hold on;
cdfplot(cell2mat(inter_all'));
cdfplot(cell2mat(intra_shuffle_all'));

%% G: anatomical cluster footprint
load([foldername_multiGeo{6},'\','neuronIndividuals_new.mat'])
[A_color1,A_color_region1]=DBSCAN_region_quantify_022422(group_ori_multiGeo{6,2},neuronIndividuals_new,[]);
load([foldername_multiGeo{5},'\','neuronIndividuals_new.mat'])
[A_color2,A_color_region2]=DBSCAN_region_quantify_022422(group_ori_multiGeo{5,2},neuronIndividuals_new,[]);

subplot(221)
imagesc(A_color1);
subplot(223)
imagesc(A_color_region1);
subplot(222)
imagesc(A_color2);
subplot(224)
imagesc(A_color_region2)

%% H: cross time stability
% output format: mice*trial cell, each trial contains 3 cells represent
% each of teh 6min period
[gp_period,A_color_period,A_color_region_period]=cross_time_stability_clust_cal(foldername_multiGeo,group_ori_multiGeo,2);

period_name={'0-6min','3-9min','6-12min'};

for tk=1:size(gp_period{6,2},2)
    subplot(3,2,2*tk-1)
    imagesc(A_color_period{6,2}{tk});
    title(period_name{tk})
    subplot(3,2,2*tk)
    imagesc(A_color_region_period{6,2}{tk});
end

%% I: region size, all 6 mice 
load('D:\Xu_clusterting_paper_prep11_2020\final_code\final_cluster_data\cluster_region_sz\fig1_3_multiGeo_cluster_Region_sz.mat')
for tk=1:6
   avg_max_region_mat(tk,:)=cell2mat(avg_region_multiGeo_k{tk,2});

   shuf_region=[];
   for i=2:10
       shuf_region(i-1)=mean(avg_region_shuf_multiGeo_k{tk,2}{i});
   end
   avg_max_reg_shuf_all_mat(tk,:)=shuf_region;
end

mean_avg_max_region_mat=mean(avg_max_region_mat,1);
se_avg_max_region_mat=std(avg_max_region_mat,[],1)/(size(avg_max_region_mat,1)^0.5);

mean_avg_max_reg_shuf_all_mat=mean(avg_max_reg_shuf_all_mat,1);
se_avg_max_reg_shuf_all_mat=std(avg_max_reg_shuf_all_mat,[],1)/(size(avg_max_reg_shuf_all_mat,1)^0.5);

x=[[0:length(mean_avg_max_region_mat)-1]',[0:length(mean_avg_max_reg_shuf_all_mat)-1]'];
y=[mean_avg_max_region_mat',mean_avg_max_reg_shuf_all_mat'];
minMaxRange={[min_avg_max_region_mat',max_avg_max_region_mat'],[min_avg_max_reg_shuf_all_mat',max_avg_max_reg_shuf_all_mat']};
% plotwithMaxMixRange(x,y,minMaxRange,{[0 0 0.9],[0 0 0.55],[0 0 0.20]},{'-o','-s','-^'});

errorbar(0:length(mean_avg_max_region_mat)-1,mean_avg_max_region_mat,se_avg_max_region_mat,'-o','color',[0    0.4470    0.7410]);
hold on;
errorbar(0:length(mean_avg_max_reg_shuf_all_mat)-1,mean_avg_max_reg_shuf_all_mat,se_avg_max_reg_shuf_all_mat,'-s','color',[0.9294    0.6941    0.1255]);

xlim([-1 9])

%% J: distance-temporal correlation plot

[all_dis_corr_circle_all_intra,all_dis_corr_circle_all_inter,all_dis_corr_circle_all_all,midx,f,gof]=pairwise_dis_tempCorr_031123(foldername_multiGeo,group_ori_multiGeo,2);

x_intra=linspace( min(all_dis_corr_circle_all_intra(:,1)), max(all_dis_corr_circle_all_intra(:,1)), 100)';
y_intra=f1.a*x_intra.^f1.b+f1.c;
x_inter=linspace( min(all_dis_corr_circle_all_inter(:,1)), max(all_dis_corr_circle_all_inter(:,1)), 100)';
y_inter=f2.a*x_inter.^f2.b+f2.c;

mean(all_dis_corr_circle_all_intra(:,1))
std(all_dis_corr_circle_all_intra(:,1))/length(all_dis_corr_circle_all_intra(:,1))^0.5

mean(all_dis_corr_circle_all_inter(:,1))
std(all_dis_corr_circle_all_inter(:,1))/length(all_dis_corr_circle_all_inter(:,1))^0.5

mean(all_dis_corr_circle_all_all(:,1))
std(all_dis_corr_circle_all_all(:,1))/length(all_dis_corr_circle_all_all(:,1))^0.5

mean(all_dis_corr_circle_all_intra(:,2))
std(all_dis_corr_circle_all_intra(:,2))/length(all_dis_corr_circle_all_intra(:,2))^0.5

mean(all_dis_corr_circle_all_inter(:,2))
std(all_dis_corr_circle_all_inter(:,2))/length(all_dis_corr_circle_all_inter(:,2))^0.5

% test the intra-cluster and inter-cluster curves
[h,p] = kstest2(y_intra,y_inter)

% spearman correlation of cyan curve 
[rho,pval] = corr(all_dis_corr_circle_all_all(1:100:end,1),all_dis_corr_circle_all_all(1:100:end,2),'Type','spearman','Rows','complete')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% main code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% finished
