run('D:\Xu_clusterting_paper_prep11_2020\final_code\data_prepare\neuron_data_info.m')
all_colors=distinguishable_colors(10);

%% A. behav
load([foldername_multiGeo{6},'\','behav.mat']);
plot(behavIndividuals{3}.position(:,1),behavIndividuals{3}.position(:,2),'color',[0.5 0.5 0.5],'lineWidth',2);

%% B: heatmap trace
load([foldername_multiGeo{6},'\','neuronIndividuals_new.mat'])
load(['D:\Xu_clusterting_paper_prep11_2020\final_code\final_cluster_data\cluster_optimal_num\Fig1_3_multiGeo_clust_original.mat'])

group_curr=group_ori_multiGeo{6,3};
clustered_timeSeries2(neuronIndividuals_new{3}.C,group_curr)

%% C: correlation map
perm=[];
for j = 1:length(unique(group_curr))
    perm = [perm;find(group_curr == j)];
end
CM22=squareform(1-pdist(neuronIndividuals_new{3}.C,'correlation'));
displayHeatmap_general(CM22,CM22,perm,all_colors,max(group_curr),group_curr,1);

%% D: Intra, inter cluster correlation, wiht shuffle, square, use cumulative distribution
for i=1:6
    load([foldername_multiGeo{i},'\','neuronIndividuals_new.mat'])
    group_curr=group_ori_multiGeo{i,3};
    [intra_all{i},inter_all{i},~,~,intra_shuffle_all{i},~]=intra_inter_cluster_corr_dis(neuronIndividuals_new,group_curr,3,'corr');
end

intra_inter_correlation_cdf(intra_all,inter_all,intra_shuffle_all)

%% E: region, footprint
%square
load([foldername_multiGeo{6},'\','neuronIndividuals_new.mat'])
[A_color1,A_color_region1]=DBSCAN_region_quantify_022422(group_ori_multiGeo{6,3},neuronIndividuals_new,[]);
load([foldername_multiGeo{3},'\','neuronIndividuals_new.mat'])
[A_color2,A_color_region2]=DBSCAN_region_quantify_022422(group_ori_multiGeo{3,3},neuronIndividuals_new,[]);

subplot(221)
imagesc(A_color1);
subplot(223)
imagesc(A_color_region1);
subplot(222)
imagesc(A_color2);
subplot(224)
imagesc(A_color_region2)

%% F: cross time stability
% output format: mice*trial cell, each trial contains 3 cells represent
% each of teh 6min period
[gp_period,A_color_period,A_color_region_period]=cross_time_stability_clust_cal(foldername_multiGeo,group_ori_multiGeo,3);

period_name={'0-6min','3-9min','6-12min'};
figure;
for tk=1:size(gp_period{6,3},2)
    subplot(3,2,2*tk-1)
    imagesc(A_color_period{6,3}{tk});
    title(period_name{tk})
    subplot(3,2,2*tk)
    imagesc(A_color_region_period{6,3}{tk});
end

%% G: cluster size across cluster num
load('D:\Xu_clusterting_paper_prep11_2020\final_code\final_cluster_data\cluster_region_sz\fig1_3_multiGeo_cluster_Region_sz.mat')
for tk=1:6
   avg_max_region_mat(tk,:)=cell2mat(avg_region_multiGeo_k{tk,3});

   shuf_region=[];
   for i=2:10
       shuf_region(i-1)=mean(avg_region_shuf_multiGeo_k{tk,3}{i});
   end
   avg_max_reg_shuf_all_mat(tk,:)=shuf_region;
end

mean_avg_max_region_mat=mean(avg_max_region_mat,1);
se_avg_max_region_mat=std(avg_max_region_mat,[],1)/(size(avg_max_region_mat,1)^0.5);

mean_avg_max_reg_shuf_all_mat=mean(avg_max_reg_shuf_all_mat,1);
se_avg_max_reg_shuf_all_mat=std(avg_max_reg_shuf_all_mat,[],1)/(size(avg_max_reg_shuf_all_mat,1)^0.5);

x=[[0:length(mean_avg_max_region_mat)-1]',[0:length(mean_avg_max_reg_shuf_all_mat)-1]'];
y=[mean_avg_max_region_mat',mean_avg_max_reg_shuf_all_mat'];
% minMaxRange={[min_avg_max_region_mat',max_avg_max_region_mat'],[min_avg_max_reg_shuf_all_mat',max_avg_max_reg_shuf_all_mat']};
% plotwithMaxMixRange(x,y,minMaxRange,{[0 0 0.9],[0 0 0.55],[0 0 0.20]},{'-o','-s','-^'});

errorbar(0:length(mean_avg_max_region_mat)-1,mean_avg_max_region_mat,se_avg_max_region_mat,'-o','color',[0    0.4470    0.7410]);
hold on;
errorbar(0:length(mean_avg_max_reg_shuf_all_mat)-1,mean_avg_max_reg_shuf_all_mat,se_avg_max_reg_shuf_all_mat,'-s','color',[0.9294    0.6941    0.1255]);

xlim([-1 9])
for i=1:9
    disp(ranksum(avg_max_region_mat(:,i),avg_max_reg_shuf_all_mat(:,i)))
end

%% H: distance-temporal correlation plot

[all_dis_corr_circle_all_intra,all_dis_corr_circle_all_inter,all_dis_corr_circle_all_all,midx,f,gof]=pairwise_dis_tempCorr_031123(foldername_multiGeo,group_ori_multiGeo,3);
f1=f{1};
f2=f{2};

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
signrank(y_intra,y_inter)

% spearman correlation of cyan curve 
[rho,pval] = corr(all_dis_corr_circle_all_all(1:100:end,1),all_dis_corr_circle_all_all(1:100:end,2),'Type','spearman','Rows','complete')


%% I: autocorrelation of cluster temporal trace
load([foldername_multiGeo{6},'\','neuronIndividuals_new.mat'])

% circle
group2=group_ori_multiGeo{6,2};

for j = 1:length(unique(group2))
    rangee=find(group2 == j);
    nC=neuronIndividuals_new{2}.C(rangee,1:end)';
    nC_sum=sum(nC,1);
    [acf{j},lag{j}] = xcorr(nC_sum','coeff');
end

for j = 1:length(unique(group2))
    subplot(length(unique(group2)),1,j);
    plot(lag{j}/15,acf{j},'-','color',all_colors(j,:));hold on; % temporal sample rate: 15Hz
    xlim([-5 5])
end

% square
group2=group_ori_multiGeo{6,3};

for j = 1:length(unique(group2))
    rangee=find(group2 == j);
    nC=neuronIndividuals_new{3}.C(rangee,1:end)';
    nC_sum=sum(nC,1);
    [acf{j},lag{j}] = xcorr(nC_sum','coeff');
end

for j = 1:length(unique(group2))
    subplot(length(unique(group2)),1,j);
    plot(lag{j}/15,acf{j},'-','color',all_colors(j,:));hold on; % temporal sample rate: 15Hz
    xlim([-5 5])
end

%% panel J-K: cluster overlap, circle and square
% circle, J
[gp_period,A_color_period,A_color_region_period]=cross_time_stability_clust_cal(foldername_multiGeo,group_ori_multiGeo,2);
overlap_circle=[];
for tk=1:6
    [~,overlap_circle(tk,1)]=new_cluster_overlap_032121(gp_period{tk,2}{1},gp_period{tk,2}{2});
    [~,overlap_circle(tk,2)]=new_cluster_overlap_032121(gp_period{tk,2}{2},gp_period{tk,2}{3});
    [~,overlap_circle(tk,3)]=new_cluster_overlap_032121(gp_period{tk,2}{1},gp_period{tk,2}{3});
end

% square, K
[gp_period,A_color_period,A_color_region_period]=cross_time_stability_clust_cal(foldername_multiGeo,group_ori_multiGeo,3);
overlap_square=[];
for tk=1:6
    [~,overlap_square(tk,1)]=new_cluster_overlap_032121(gp_period{tk,3}{1},gp_period{tk,3}{2});
    [~,overlap_square(tk,2)]=new_cluster_overlap_032121(gp_period{tk,3}{2},gp_period{tk,3}{3});
    [~,overlap_square(tk,3)]=new_cluster_overlap_032121(gp_period{tk,3}{1},gp_period{tk,3}{3});
end

%% panel L,M: avg dis inside 35um/outside 35um
% circle
[all_dis_corr_circle_all_intra,all_dis_corr_circle_all_inter]=pairwise_dis_tempCorr_031123(foldername_multiGeo,group_ori_multiGeo,2);

mice_circle_dis_corr_intra_close=[];
mice_circle_dis_corr_intra_far=[]; 
mice_circle_dis_corr_inter_close=[];
mice_circle_dis_corr_inter_far=[]; 
for tk=1:6
    mice_circle_dis_corr_intra=all_dis_corr_circle_all_intra(midx{1}==tk,:);
    mice_circle_dis_corr_intra_close(tk,1)=mean(mice_circle_dis_corr_intra(mice_circle_dis_corr_intra(:,1)<=35,2)); % we reverse the spatial downsample in pairwise_dis_tempCorr, so now the thres in 35
    mice_circle_dis_corr_intra_far(tk,1)=mean(mice_circle_dis_corr_intra(mice_circle_dis_corr_intra(:,1)>35,2));
    mice_circle_dis_corr_inter=all_dis_corr_circle_all_inter(midx{2}==tk,:);
    mice_circle_dis_corr_inter_close(tk,1)=mean(mice_circle_dis_corr_inter(mice_circle_dis_corr_inter(:,1)<=35,2));
    mice_circle_dis_corr_inter_far(tk,1)=mean(mice_circle_dis_corr_inter(mice_circle_dis_corr_inter(:,1)>35,2));
end

%square
[all_dis_corr_square_all_intra,all_dis_corr_square_all_inter]=pairwise_dis_tempCorr_031123(foldername_multiGeo,group_ori_multiGeo,3);

mice_square_dis_corr_intra_close=[];
mice_square_dis_corr_intra_far=[]; 
mice_square_dis_corr_inter_close=[];
mice_square_dis_corr_inter_far=[]; 
for tk=1:6
    mice_square_dis_corr_intra=all_dis_corr_square_all_intra(midx_square{1}==tk,:);
    mice_square_dis_corr_intra_close(tk,1)=mean(mice_square_dis_corr_intra(mice_square_dis_corr_intra(:,1)<=35,2));
    mice_square_dis_corr_intra_far(tk,1)=mean(mice_square_dis_corr_intra(mice_square_dis_corr_intra(:,1)>35,2));
    mice_square_dis_corr_inter=all_dis_corr_square_all_inter(midx_square{2}==tk,:);
    mice_square_dis_corr_inter_close(tk,1)=mean(mice_square_dis_corr_inter(mice_square_dis_corr_inter(:,1)<=35,2));
    mice_square_dis_corr_inter_far(tk,1)=mean(mice_square_dis_corr_inter(mice_square_dis_corr_inter(:,1)>35,2));
end

%% N,O: cluster size, optimal num
% circle
avg_region=[];
avg_region_shuf=[];
for i=1:6
    load([foldername_multiGeo{i},'\','neuronIndividuals_new.mat'])
    [~,~,avg_region(i,1),~,~,~,avg_region_shuf(i,1)]=DBSCAN_region_quantify_022422(group_ori_multiGeo{i,2},neuronIndividuals_new,[]);
end
ranksum(avg_region,avg_region_shuf)

% square 
avg_region=[];
avg_region_shuf=[];
for i=1:6
    load([foldername_multiGeo{i},'\','neuronIndividuals_new.mat'])
    [~,~,avg_region(i,1),~,~,~,avg_region_shuf(i,1)]=DBSCAN_region_quantify_022422(group_ori_multiGeo{i,3},neuronIndividuals_new,[]);
end
ranksum(avg_region,avg_region_shuf)
