%% Supplemental Fig6: supporting data for Fig2 and suppl5

run('D:\Xu_clusterting_paper_prep11_2020\final_code\data_prepare\neuron_data_info.m')
all_colors=distinguishable_colors(10);
load(['D:\Xu_clusterting_paper_prep11_2020\final_code\final_cluster_data\cluster_optimal_num\Fig2_cir_rec_clust_original.mat'])

foldername=foldername_fig2;
%% panel A,E: behav trajectory
% circle
load([foldername{1},'\','Behav.mat']);
figure;
plot(behavIndividuals{2}.position(:,1),behavIndividuals{2}.position(:,2),'color',[0.5,0.5,0.5],'lineWidth',2);

% square
figure;
plot(behavIndividuals{1}.position(:,1),behavIndividuals{1}.position(:,2),'color',[0.5,0.5,0.5],'lineWidth',2);


%% panel B,F: intra, intra shuffle and inter clust corr, all mice
% circle
for tk=1:12
    cond_sign=1;
    load([foldername{tk},'\','neuronIndividuals_new.mat'])
    group_model=group_ori_fig2{tk,cond_sign};
    [intra_all{tk},inter_all{tk},~,~,intra_shuffle_all{tk},~]=intra_inter_cluster_corr_dis(neuronIndividuals_new,group_model,cond_sign,'corr');
end

colorClusters_all=colormap(lines);close

a1x=figure;

h1=cdfplot(cell2mat(intra_all'));
hold on;
h2=cdfplot(cell2mat(inter_all'));
h3=cdfplot(cell2mat(intra_shuffle_all'));
set(h1,'color',[0.4940,0.1840,0.5560]);
set(h2,'color',[0.4660,0.6740,0.1880]);
set(h3,'color',[0.3010,0.7450,0.9330]);

p1=infer_cdf_loc(cell2mat(intra_all'),mean(cell2mat(intra_all')));
p2=infer_cdf_loc(cell2mat(inter_all'),mean(cell2mat(inter_all')));
p3=infer_cdf_loc(cell2mat(intra_shuffle_all'),mean(cell2mat(intra_shuffle_all')));

plot(p1(1),p1(2),'.','color',colorClusters_all(4,:),'MarkerSize',20)
plot(p2(1),p2(2),'.','color',colorClusters_all(5,:),'MarkerSize',20)
plot(p3(1),p3(2),'.','color',colorClusters_all(6,:),'MarkerSize',20)

% square
for tk=1:12
    cond_sign=2;
    load([foldername{tk},'\','neuronIndividuals_new.mat'])
    group_model=group_ori_fig2{tk,cond_sign};
    [intra_all{tk},inter_all{tk},~,~,intra_shuffle_all{tk},~]=intra_inter_cluster_corr_dis(neuronIndividuals_new,group_model,cond_sign,'corr');
end

colorClusters_all=colormap(lines);close

a1x=figure;

h1=cdfplot(cell2mat(intra_all'));
hold on;
h2=cdfplot(cell2mat(inter_all'));
h3=cdfplot(cell2mat(intra_shuffle_all'));
set(h1,'color',[0.4940,0.1840,0.5560]);
set(h2,'color',[0.4660,0.6740,0.1880]);
set(h3,'color',[0.3010,0.7450,0.9330]);

p1=infer_cdf_loc(cell2mat(intra_all'),mean(cell2mat(intra_all')));
p2=infer_cdf_loc(cell2mat(inter_all'),mean(cell2mat(inter_all')));
p3=infer_cdf_loc(cell2mat(intra_shuffle_all'),mean(cell2mat(intra_shuffle_all')));

plot(p1(1),p1(2),'.','color',colorClusters_all(4,:),'MarkerSize',20)
plot(p2(1),p2(2),'.','color',colorClusters_all(5,:),'MarkerSize',20)
plot(p3(1),p3(2),'.','color',colorClusters_all(6,:),'MarkerSize',20)

%% panel C,G: distance-temporal correlation plot, all mice

% circle
[all_dis_corr_circle_all_intra,all_dis_corr_circle_all_inter,all_dis_corr_circle_all_all,midx,f,gof]=pairwise_dis_tempCorr_031123(foldername_fig2,group_ori_fig2,1);
[rho,pval] = corr(all_dis_corr_circle_all_all,all_dis_corr_circle_all_all,'Type','spearman','Rows','complete');

% square
[all_dis_corr_square_all_intra,all_dis_corr_square_all_inter,all_dis_corr_square_all_all,midx,f,gof]=pairwise_dis_tempCorr_031123(foldername_fig2,group_ori_fig2,1);
[rho,pval] = corr(all_dis_corr_square_all_all,all_dis_corr_square_all_all,'Type','spearman','Rows','complete');

%% panel D,H: cluster region size change with cluster number, all mice
load('D:\Xu_clusterting_paper_prep11_2020\final_code\final_cluster_data\cluster_region_sz\fig2_cir_rec_cluster_Region_sz.mat')

% circle
avg_max_region_mat=[];
avg_max_reg_shuf_all_mat=[];
for tk=1:6
   avg_max_region_mat(tk,:)=cell2mat(avg_region_fig2_k{tk,1});
   shuf_region=[];
   for i=2:10
       shuf_region(i-1)=mean(avg_region_shuf_fig2_k{tk,1}{i});
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

% square
avg_max_region_mat=[];
avg_max_reg_shuf_all_mat=[];
for tk=1:6
   avg_max_region_mat(tk,:)=cell2mat(avg_region_fig2_k{tk,2});
   shuf_region=[];
   for i=2:10
       shuf_region(i-1)=mean(avg_region_shuf_fig2_k{tk,2}{i});
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

