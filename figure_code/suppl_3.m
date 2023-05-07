run('D:\Xu_clusterting_paper_prep11_2020\final code\data_prepare\data_info.m')
all_colors=distinguishable_colors(10);

%% A. behav
load([foldername_multiGeo{6},'\','square1_behav.mat']);
plot(behav.position(:,1),behav.position(:,2),'color',[0.5 0.5 0.5],'lineWidth',2);

%% B: heatmap trace
load([foldername_multiGeo{6},'\','neuronIndividuals_new.mat'])
load(['D:\Xu_clusterting_paper_prep11_2020\final data\final_cluster_data\Fig1_3_multiGeo_clust_original.mat'])

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

%% C: autocorrelation of cluster temporal trace
for j = 1:length(unique(group2))
    rangee=find(group2 == j);
    nC=neuronIndividuals_new{3}.C(rangee,1:end)';
    nC_sum=sum(nC,1);
    [acf{j},lag{j}] = xcorr(nC_sum','coeff');
end

for j = 1:length(unique(group2))
    subplot(length(unique(group2)),1,j);
    plot(lag{j},acf{j},'-','color',colorClusters2(j,:));hold on;
    xlim([-75 75])
    ctt=ctt+1;
end


%% Panel d: region, footprint
%square
load('C:\Users\exx\Desktop\HDAC paper fig and method\SFN2019\SFN2019 fig\SCI clustNum plot\tri cir sqr\tri_cir_sqr_dataset.mat');

load(['square1_placement_cell_cluster','\','variables_clustering_0.35PC.mat'])
group_model=group{3};

% DBSCAN_region_quantify_func_simplify(gp_rec,clust_num,shuffles_num,clust_idx-1,cond_sign,neuronIndividuals_new)
[A_color1,A_color_region1,avg_max_region,avg_max_reg_shuf_all,avg_max_reg_shuf_all_95,A_color_rand,A_color_region_rand]=DBSCAN_region_quantify_func_single_group(group_model,3,neuronIndividuals_new)

load('D:\Remapping_square_circle_triangle_061119_061319\M3421F\neuronIndividuals_new.mat')
load(['D:\Remapping_square_circle_triangle_061119_061319\M3421F\square1_placement_cell_cluster','\','variables_clustering_0.35PC.mat'])
group_model=group{3};

[A_color2,A_color_region2,avg_max_region,avg_max_reg_shuf_all,avg_max_reg_shuf_all_95,A_color_rand,A_color_region_rand]=DBSCAN_region_quantify_func_single_group(group_model,3,neuronIndividuals_new)

subplot(221)
imagesc(A_color1);
subplot(223)
imagesc(A_color_region1);
subplot(222)
imagesc(A_color2);
subplot(224)
imagesc(A_color_region2)
%% Panel e: Intra, inter cluster correlation, wiht shuffle, circle
for tk=1:6
    cond_sign=3;
    load([foldername{tk},'\','neuronIndividuals_new.mat'])
    load([foldername{tk},'\','square1_placement_cell_cluster','\','variables_clustering_0.35PC.mat'])
    group_model=group{3};
    [intra_all{tk},inter_all{tk},~,~,intra_shuffle_all{tk},~]=intra_inter_cluster_corr_dis(neuronIndividuals_new,group_model,cond_sign,'corr');
end
colorClusters_all=colormap(lines);close

a1x=figure;
% subplot(211)

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

ranksum(cell2mat(inter_all'),cell2mat(intra_shuffle_all'))

mean(cell2mat(intra_all'))
std(cell2mat(intra_all'))/length(cell2mat(intra_all'))^0.5
mean(cell2mat(inter_all'))
std(cell2mat(inter_all'))/length(cell2mat(inter_all'))^0.5
mean(cell2mat(intra_shuffle_all'))
std(cell2mat(intra_shuffle_all'))/length(cell2mat(intra_shuffle_all'))^0.5

% subplot(212)
% b1=bar([mean(cell2mat(intra_all')),mean(cell2mat(inter_all')),mean(cell2mat(intra_shuffle_all'));0 0 0]);
% for i=1:3
%     b1(i).FaceColor=[colorClusters_all(i+3,:)]; %in 2017a
% %     b1.CData(i,:)=[colorClusters_all(i,:)] %in 2017b and above
% end
k1=cell2mat(intra_all');
k2=cell2mat(inter_all');
k3=cell2mat(intra_shuffle_all');

[H1,P1] = kstest2(k1(1:20:end),k2(1:100:end));
[H2,P2] = kstest2(k1(1:20:end),k3(1:10000:end));
[H3,P3] = kstest2(k2(1:100:end),k3(1:10000:end));


%% panel f: distance-temporal correlation plot
mice_all=[1:6];
cluster_file_name={};
for tk=1:6
    fname=foldername{mice_all(tk)};
    cluster_file_name{tk}=[fname,'\','square1_placement_cell_cluster','\','variables_clustering_0.35PC.mat'];
    session=3;
end
% [all_dis_corr_circle_all_intra,all_dis_corr_circle_all_inter,~,f1,f2]=pairwise_dis_tempCorr(foldername,cluster_file_name,session);
[all_dis_corr_circle_all_intra,all_dis_corr_circle_all_inter,all_dis_corr_circle_all_all,midx,f,gof]=pairwise_dis_tempCorr_091421(foldername,cluster_file_name,session);

[H, pValue, KSstatistic] = kstest_2s_2d(all_dis_corr_circle_all_intra(1:1000:end,:), all_dis_corr_circle_all_inter(1:1000:end,:));

f1=f{1};
f2=f{2};
x_intra=linspace( min(all_dis_corr_circle_all_intra(:,1)), max(all_dis_corr_circle_all_intra(:,1)), 100)';
y_intra=f1.a*x_intra.^f1.b+f1.c;
x_inter=linspace( min(all_dis_corr_circle_all_inter(:,1)), max(all_dis_corr_circle_all_inter(:,1)), 100)';
y_inter=f2.a*x_inter.^f2.b+f2.c;

mean(y_intra)
std(y_intra)/length(y_intra)^0.5
mean(y_inter)
std(y_inter)/length(y_inter)^0.5

mean(all_dis_corr_circle_all_intra(:,1)) % distance here already been corrected
std(all_dis_corr_circle_all_intra(:,1))/length(all_dis_corr_circle_all_intra(:,1))^0.5

mean(all_dis_corr_circle_all_inter(:,1))
std(all_dis_corr_circle_all_inter(:,1))/length(all_dis_corr_circle_all_inter(:,1))^0.5

mean(all_dis_corr_circle_all_all(:,1))
std(all_dis_corr_circle_all_all(:,1))/length(all_dis_corr_circle_all_all(:,1))^0.5

mean(all_dis_corr_circle_all_intra(:,2))
std(all_dis_corr_circle_all_intra(:,2))/length(all_dis_corr_circle_all_intra(:,2))^0.5

mean(all_dis_corr_circle_all_inter(:,2))
std(all_dis_corr_circle_all_inter(:,2))/length(all_dis_corr_circle_all_inter(:,2))^0.5

[h,p] = kstest2(y_intra,y_inter)
ranksum(y_intra,y_inter)

[rho,pval] = corr(all_dis_corr_circle_all_all(1:100:end,1),all_dis_corr_circle_all_all(1:100:end,2),'Type','spearman','Rows','complete')

%% panel g: region size, all 6 mice (RESULTS ARE SAVED FROM ROUND 10). CAN DIRECTLY BE LOADED
avg_max_region={};
avg_max_reg_shuf_all={};
avg_max_reg_shuf_all_95={};
for tk=1:6
    load([foldername{tk},'\neuronIndividuals_new.mat']);
    gp_rec=group_record(:,tk)
    clust_num=[2:10];
    shuffles_num=100;
    cond_sign=3;

    load([foldername{tk},'\square1_placement_cell_cluster','\','variables_clustering_0.35PC.mat'])
    group_model=group{3};
    clust_idx=length(unique(group_model));
    gp_rec{cond_sign}{clust_idx-1}=group_model;

    [~,~,avg_max_region{tk},avg_max_reg_shuf_all{tk},avg_max_reg_shuf_all_95{tk}]=DBSCAN_region_quantify_func_simplify(gp_rec,clust_num,shuffles_num,clust_idx-1,cond_sign,neuronIndividuals_new);
    close
end

for tk=1:6
   avg_max_region_mat(tk,:)=avg_max_region{tk}(3,:);
   avg_max_reg_shuf_all_mat(tk,:)=avg_max_reg_shuf_all{tk}(3,:);
   avg_max_reg_shuf_all_95_mat(tk,:)=avg_max_reg_shuf_all_95{tk}(3,:);
end

mean_avg_max_region_mat=mean(avg_max_region_mat,1);
se_avg_max_region_mat=std(avg_max_region_mat,[],1)/(size(avg_max_region_mat,1)^0.5);
min_avg_max_region_mat=mean(avg_max_region_mat,1)-se_avg_max_region_mat;
max_avg_max_region_mat=mean(avg_max_region_mat,1)+se_avg_max_region_mat;

mean_avg_max_reg_shuf_all_mat=mean(avg_max_reg_shuf_all_mat,1);
se_avg_max_reg_shuf_all_mat=std(avg_max_reg_shuf_all_mat,[],1)/(size(avg_max_reg_shuf_all_mat,1)^0.5);
min_avg_max_reg_shuf_all_mat=mean(avg_max_reg_shuf_all_mat,1)-se_avg_max_reg_shuf_all_mat;
max_avg_max_reg_shuf_all_mat=mean(avg_max_reg_shuf_all_mat,1)+se_avg_max_reg_shuf_all_mat;

mean_avg_max_reg_shuf_all_95_mat=mean(avg_max_reg_shuf_all_95_mat,1);
se_avg_max_reg_shuf_all_95_mat=std(avg_max_reg_shuf_all_95_mat,[],1)/(size(avg_max_reg_shuf_all_95_mat,1)^0.5);
min_avg_max_reg_shuf_all_95_mat=mean(avg_max_reg_shuf_all_95_mat,1)-se_avg_max_reg_shuf_all_95_mat;
max_avg_max_reg_shuf_all_95_mat=mean(avg_max_reg_shuf_all_95_mat,1)+se_avg_max_reg_shuf_all_95_mat;

x=[[0:length(mean_avg_max_region_mat)-1]',[0:length(mean_avg_max_reg_shuf_all_mat)-1]',[0:length(mean_avg_max_reg_shuf_all_95_mat)-1]'];
y=[mean_avg_max_region_mat',mean_avg_max_reg_shuf_all_mat',mean_avg_max_reg_shuf_all_95_mat'];
minMaxRange={[min_avg_max_region_mat',max_avg_max_region_mat'],[min_avg_max_reg_shuf_all_mat',max_avg_max_reg_shuf_all_mat'],[min_avg_max_reg_shuf_all_95_mat',max_avg_max_reg_shuf_all_95_mat']};
% plotwithMaxMixRange(x,y,minMaxRange,{[0 0 0.9],[0 0 0.55],[0 0 0.20]},{'-o','-s','-^'});

errorbar(0:length(mean_avg_max_region_mat)-1,mean_avg_max_region_mat,se_avg_max_region_mat,'-o','color',[0    0.4470    0.7410]);
hold on;
errorbar(0:length(mean_avg_max_reg_shuf_all_mat)-1,mean_avg_max_reg_shuf_all_mat,se_avg_max_reg_shuf_all_mat,'-s','color',[0.8500    0.3250    0.0980]);
errorbar(0:length(mean_avg_max_reg_shuf_all_95_mat)-1,mean_avg_max_reg_shuf_all_95_mat,se_avg_max_reg_shuf_all_95_mat,'-^','color',[0.9290    0.6940    0.1250]);

xlim([-1 9])

%% panel h cross time stability
group_all_mice={};
adjrand_idx_val_mice={};
adjrand_idx_val_shuffle_mice={};

for tk=1:length(foldername)
    load([foldername{tk},'\','neuronIndividuals_new.mat']);
    win_leng=5400;%4min
    session=3;
    time_series=1:win_leng/2:size(neuronIndividuals_new{session}.C,2)-win_leng;

    load([foldername{tk},'\','square1_placement_cell_cluster','\','variables_clustering_0.35PC.mat'])
    group_model=group{3}; 

    [~,group_all_all{tk}]=cluster_stability_illustrate_func(time_series,session,neuronIndividuals_new,win_leng,group_model);
end
for tk=1:length(foldername)
    load([foldername{tk},'\','neuronIndividuals_new.mat']);
    group_all=group_all_all{tk};
    figure;
    for i=1:length(group_all)
        subplot(2,9,i);
        [A_color,A_color_region,avg_max_region,avg_max_reg_shuf_all,avg_max_reg_shuf_all_95]=DBSCAN_region_quantify_func_single_group_no_plot(group_all{i},neuronIndividuals_new);
        imagesc(A_color);
        subplot(2,9,i+9);
        imagesc(A_color_region);
    end
end


adjoverlap_val_mice={};
overlap_val_mice={};
adjoverlap_val_shuffle_mice={};
overlap_val_shuffle_mice={};
adjoverlap_val_shuffle_95_mice={};
overlap_val_shuffle_95_mice={};
overlap_val_shuffle_nc_mice={};

for tk=1:length(foldername)
    group_all=group_all_all{tk};
    adjoverlap_val=[];
    overlap_val=[];
    nc=[];
    
    adjoverlap_shuffle_val=[];
    overlap_shuffle_val=[];
    
    for i=1:length(group_all)-1
%         [adjoverlap_val(i),overlap_val(i),~,~,nc(i)]=RandIndex(group_all{i},group_all{i+1}); % this calculation is correct, proved by sklearn ari calculation
        overlap_val(i)=new_cluster_overlap(group_all{i},group_all{i+1});
    end
%     [adjoverlap_val(i+1),overlap_val(i+1),~,~,nc(i+1)]=RandIndex(group_all{3},group_all{1});
    overlap_val(i+1)=new_cluster_overlap(group_all{3},group_all{1});
    
    shuffle=1000;
    for i=1:length(group_all)-1
        for j=1:shuffle
            g1=group_all{i}(randperm(length(group_all{i})));
            g2=group_all{i+1}(randperm(length(group_all{i+1})));
%             [adjoverlap_shuffle_val(i,j),overlap_shuffle_val(i,j)]=RandIndex(g1,g2);
            overlap_shuffle_val(i,j)=new_cluster_overlap(g1,g2);
        end
    end
    
    for j=1:shuffle
        g3=group_all{3}(randperm(length(group_all{3})));
        g4=group_all{1}(randperm(length(group_all{1})));
%         [adjoverlap_shuffle_val(i+1,j),overlap_shuffle_val(i+1,j)]=RandIndex(g1,g2);
        overlap_shuffle_val(i+1,j)=new_cluster_overlap(g1,g2);
    end
    
    adjoverlap_val_mice{tk}=adjoverlap_val;
    overlap_val_mice{tk}=overlap_val;
    adjoverlap_shuffle_val_mice{tk}=adjoverlap_shuffle_val;
    overlap_shuffle_val_mice{tk}=overlap_shuffle_val;   
    adjoverlap_val_shuffle_qt_mice{tk}=quantile(adjoverlap_shuffle_val,.95,2);
    overlap_val_shuffle_qt_mice{tk}=quantile(overlap_shuffle_val,.95,2);
    overlap_val_shuffle_nc_mice{tk}=nc;
%     adjoverlap_val_shuffle_mice{tk}=adjoverlap_val_shuffle;
end


ori_adjrand=cell2mat(adjoverlap_val_mice');
ori_rand=cell2mat(overlap_val_mice');
ori_rand_nc=cell2mat(overlap_val_shuffle_nc_mice');

ori_adjrand_shuffle_qt=cell2mat(adjoverlap_val_shuffle_qt_mice)';
ori_rand_shuffle_qt=cell2mat(overlap_val_shuffle_qt_mice)';


%% panel j new: cluster size--M3425, M3421
load('C:\Users\exx\Desktop\HDAC paper fig and method\SFN2019\SFN2019 fig\SCI clustNum plot\tri cir sqr\tri_cir_sqr_dataset.mat');

mice_all=[1:6];
cluster_file_name={};
avg_max_region={};
avg_max_reg_shuf_all={};
avg_max_reg_shuf_all_95={};
for tk=[6 3]
    figure;
    load([foldername{tk},'\neuronIndividuals_new.mat']);
    gp_rec=group_record(:,tk);
    clust_num=[2:10];
    shuffles_num=100;
    cond_sign=3;

    load([foldername{tk},'\square1_placement_cell_cluster','\','variables_clustering_0.35PC.mat'])
    group_model=group{3};
    clust_idx=length(unique(group_model));
    gp_rec{cond_sign}{clust_idx-1}=group_model;

    [~,~,avg_max_region{tk},avg_max_reg_shuf_all{tk},avg_max_reg_shuf_all_95{tk}]=DBSCAN_region_quantify_func_simplify_no_plot(gp_rec,clust_num,shuffles_num,clust_idx-1,cond_sign,neuronIndividuals_new);


    avg_max_region_mat=avg_max_region{tk}(cond_sign,:);
    avg_max_reg_shuf_all_mat=avg_max_reg_shuf_all{tk}(cond_sign,:);
    avg_max_reg_shuf_all_95_mat=avg_max_reg_shuf_all_95{tk}(cond_sign,:);

    mean_avg_max_region_mat=mean(avg_max_region_mat,1);
    se_avg_max_region_mat=std(avg_max_region_mat,[],1)/(size(avg_max_region_mat,1)^0.5);
    min_avg_max_region_mat=mean(avg_max_region_mat,1)-se_avg_max_region_mat;
    max_avg_max_region_mat=mean(avg_max_region_mat,1)+se_avg_max_region_mat;

    mean_avg_max_reg_shuf_all_mat=mean(avg_max_reg_shuf_all_mat,1);
    se_avg_max_reg_shuf_all_mat=std(avg_max_reg_shuf_all_mat,[],1)/(size(avg_max_reg_shuf_all_mat,1)^0.5);
    min_avg_max_reg_shuf_all_mat=mean(avg_max_reg_shuf_all_mat,1)-se_avg_max_reg_shuf_all_mat;
    max_avg_max_reg_shuf_all_mat=mean(avg_max_reg_shuf_all_mat,1)+se_avg_max_reg_shuf_all_mat;

    mean_avg_max_reg_shuf_all_95_mat=mean(avg_max_reg_shuf_all_95_mat,1);
    se_avg_max_reg_shuf_all_95_mat=std(avg_max_reg_shuf_all_95_mat,[],1)/(size(avg_max_reg_shuf_all_95_mat,1)^0.5);
    min_avg_max_reg_shuf_all_95_mat=mean(avg_max_reg_shuf_all_95_mat,1)-se_avg_max_reg_shuf_all_95_mat;
    max_avg_max_reg_shuf_all_95_mat=mean(avg_max_reg_shuf_all_95_mat,1)+se_avg_max_reg_shuf_all_95_mat;

    errorbar(0:length(mean_avg_max_region_mat)-1,mean_avg_max_region_mat,se_avg_max_region_mat,'-o','color',[0    0.4470    0.7410]);
    hold on;
    errorbar(0:length(mean_avg_max_reg_shuf_all_mat)-1,mean_avg_max_reg_shuf_all_mat,se_avg_max_reg_shuf_all_mat,'-s','color',[0.8500    0.3250    0.0980]);
    errorbar(0:length(mean_avg_max_reg_shuf_all_95_mat)-1,mean_avg_max_reg_shuf_all_95_mat,se_avg_max_reg_shuf_all_95_mat,'-^','color',[0.9290    0.6940    0.1250]);

    xlim([-1 9])
    ylim([0 14000]);
end

%% panel k: avg dis inside 35um/outside 35um
% circle
mice_all=[1:6];
cluster_file_name={};
for tk=1:6
    fname=foldername{mice_all(tk)};
    cluster_file_name{tk}=[fname,'\','circle1_placement_cell_cluster','\','variables_clustering_0.35PC.mat'];
    session=2;
end
[all_dis_corr_circle_all_intra,all_dis_corr_circle_all_inter,midx]=pairwise_dis_tempCorr(foldername,cluster_file_name,session);
 
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
mice_all=[1:6];
cluster_file_name={};
for tk=1:6
    fname=foldername{mice_all(tk)};
    cluster_file_name{tk}=[fname,'\','square1_placement_cell_cluster','\','variables_clustering_0.35PC.mat'];
    session=3;
end
[all_dis_corr_square_all_intra,all_dis_corr_square_all_inter,midx_square]=pairwise_dis_tempCorr(foldername,cluster_file_name,session);
 
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

%% panel m and n: half length of auto correlation, average dis in hlac, circle/square
%auto correlation
acf={};
 lag={};
 session=2;
 for i=1:length(foldername)
     load([foldername{i},'\','neuronIndividuals_new.mat'])
     load([foldername{i},'\','circle1_placement_cell_cluster','\','variables_clustering_0.35PC.mat'])
     group_model=group{ session};
     
    for j = 1:length(unique(group_model))
        rangee=find(group_model == j);
        nC=neuronIndividuals_new{session}.C(rangee,:)';
        nC_sum=sum(nC,1);
        [acf{i}{j},lag{i}{j}] = xcorr(nC_sum','coeff');
    end
 end
 
 m_hlhh_circle=[];
 for i=1:length(acf)
     acft=acf{i};
     acft_hlhh=[];
     for j=1:length(acft)
         acft_t=acft{j};
         max_acft_t=max(acft_t);
%          acft_t_tlength=length(find(acft_t>0.5*max_acft_t));
         acft_t_tlength=length(acft_t)/2;
         acft_hlhh(j)=acft_t_tlength;
     end
     m_hlhh_circle(i)=mean(acft_hlhh);
 end
 
 m_hlhh_circle=(m_hlhh_circle/15)';
 
 acf={};
 lag={};
 session=3;
 for i=1:length(foldername)
     load([foldername{i},'\','neuronIndividuals_new.mat'])
     load([foldername{i},'\','square1_placement_cell_cluster','\','variables_clustering_0.35PC.mat'])
     group_model=group{ session};
     
    for j = 1:length(unique(group_model))
        rangee=find(group_model == j);
        nC=neuronIndividuals_new{session}.C(rangee,:)';
        nC_sum=sum(nC,1);
        [acf{i}{j},lag{i}{j}] = xcorr(nC_sum','coeff');
    end
 end
 
 m_hlhh_square=[];
 for i=1:length(acf)
     acft=acf{i};
     acft_hlhh=[];
     for j=1:length(acft)
         acft_t=acft{j};
         max_acft_t=max(acft_t);
%          acft_t_tlength=length(find(acft_t>0.5*max_acft_t));
         acft_t_tlength=length(acft_t)/2;
         acft_hlhh(j)=acft_t_tlength;
     end
     m_hlhh_square(i)=mean(acft_hlhh);
 end
 
 m_hlhh_square=(m_hlhh_square/15)';
 
distance_sec_cir_mean=[];
distance_sec_sqr_mean=[];
cir_radius=zeros(1,length(foldername));
sqr_radius=zeros(1,length(foldername));
for i=1:length(foldername)
    t_leng_c=m_hlhh_circle(i)*15;
    load([foldername{i},'\','circle1_Behav.mat']);
    dis_all_sec=[];
    for j=1:size(behav.position,1)-t_leng_c
        position_sec=behav.position(j:j+t_leng_c-1,:);
        distance_sec=0;
        for k=1:size(position_sec,1)-1
            distance_sec=distance_sec+nansum((position_sec(k+1,:)-position_sec(k,:)).^2).^0.5;
        end
        dis_sec_all(j)=distance_sec;
    end
    distance_sec_cir_mean(i)=nanmean(dis_sec_all);
    
    cir_radius(i)=mean([max(behav.position(:,1))-min(behav.position(:,1)),max(behav.position(:,2))-min(behav.position(:,2))]);

    
    load([foldername{i},'\','square1_Behav.mat']);
    t_leng_s=m_hlhh_square(i)*15;
    dis_sec_all=[];
    for j=1:size(behav.position,1)-t_leng_s
        position_sec=behav.position(j:j+t_leng_s-1,:);
        distance_sec=0;
        for k=1:size(position_sec,1)-1
            distance_sec=distance_sec+nansum((position_sec(k+1,:)-position_sec(k,:)).^2).^0.5;
        end
        dis_sec_all(j)=distance_sec;
    end
    distance_sec_sqr_mean(i)=nanmean(dis_sec_all);
    sqr_radius(i)=mean([max(behav.position(:,1))-min(behav.position(:,1)),max(behav.position(:,2))-min(behav.position(:,2))]);

end

distance_sec_sqr_mean=distance_sec_sqr_mean';
distance_sec_cir_mean=distance_sec_cir_mean';

%suppl 1 new

%% panel k: statistics of pairwise corr for very close vs other distance range, all mice per mice
% circle
mice_all=[1:6];
cluster_file_name={};
for tk=1:6
    fname=foldername{mice_all(tk)};
    cluster_file_name{tk}=[fname,'\','circle1_placement_cell_cluster','\','variables_clustering_0.35PC.mat'];
    session=2;
end
[all_dis_corr_circle_all_intra,all_dis_corr_circle_all_inter,midx]=pairwise_dis_tempCorr(foldername,cluster_file_name,session);
 
mice_circle_dis_corr_intra_close=[];
mice_circle_dis_corr_intra_far=[]; 
mice_circle_dis_corr_inter_close=[];
mice_circle_dis_corr_inter_far=[]; 
for tk=1:6
    mice_circle_dis_corr_intra=all_dis_corr_circle_all_intra(midx{1}==tk,:);
    mice_circle_dis_corr_intra_close(tk,1)=mean(mice_circle_dis_corr_intra(mice_circle_dis_corr_intra(:,1)<=25,2));
    mice_circle_dis_corr_intra_far(tk,1)=mean(mice_circle_dis_corr_intra(mice_circle_dis_corr_intra(:,1)>25,2));
    mice_circle_dis_corr_inter=all_dis_corr_circle_all_inter(midx{2}==tk,:);
    mice_circle_dis_corr_inter_close(tk,1)=mean(mice_circle_dis_corr_inter(mice_circle_dis_corr_inter(:,1)<=25,2));
    mice_circle_dis_corr_inter_far(tk,1)=mean(mice_circle_dis_corr_inter(mice_circle_dis_corr_inter(:,1)>25,2));
end

%square
mice_all=[1:6];
cluster_file_name={};
for tk=1:6
    fname=foldername{mice_all(tk)};
    cluster_file_name{tk}=[fname,'\','square1_placement_cell_cluster','\','variables_clustering_0.35PC.mat'];
    session=3;
end
[all_dis_corr_square_all_intra,all_dis_corr_square_all_inter,midx_square]=pairwise_dis_tempCorr(foldername,cluster_file_name,session);
 
mice_square_dis_corr_intra_close=[];
mice_square_dis_corr_intra_far=[]; 
mice_square_dis_corr_inter_close=[];
mice_square_dis_corr_inter_far=[]; 
for tk=1:6
    mice_square_dis_corr_intra=all_dis_corr_square_all_intra(midx_square{1}==tk,:);
    mice_square_dis_corr_intra_close(tk,1)=mean(mice_square_dis_corr_intra(mice_square_dis_corr_intra(:,1)<=25,2));
    mice_square_dis_corr_intra_far(tk,1)=mean(mice_square_dis_corr_intra(mice_square_dis_corr_intra(:,1)>25,2));
    mice_square_dis_corr_inter=all_dis_corr_square_all_inter(midx_square{2}==tk,:);
    mice_square_dis_corr_inter_close(tk,1)=mean(mice_square_dis_corr_inter(mice_square_dis_corr_inter(:,1)<=25,2));
    mice_square_dis_corr_inter_far(tk,1)=mean(mice_square_dis_corr_inter(mice_square_dis_corr_inter(:,1)>25,2));
end

mean(mice_circle_dis_corr_intra_close)

ranksum(mice_circle_dis_corr_intra_close,mice_circle_dis_corr_intra_far)
ranksum(mice_circle_dis_corr_inter_close,mice_circle_dis_corr_inter_far)
ranksum(mice_square_dis_corr_intra_close,mice_circle_dis_corr_intra_far)
ranksum(mice_square_dis_corr_inter_close,mice_circle_dis_corr_inter_far)
ranksum(mice_square_dis_corr_intra_far,mice_circle_dis_corr_inter_far)
ranksum(mice_square_dis_corr_intra_close,mice_circle_dis_corr_inter_close)

%% panel l: original rand index, circle square
% circle
for tk=1:length(foldername)
    load([foldername{tk},'\','neuronIndividuals_new.mat']);
    win_leng=5400;%6min
    session=2;
    time_series=1:win_leng/2:size(neuronIndividuals_new{session}.C,2)-win_leng;
    load([foldername{tk},'\','circle1_placement_cell_cluster','\','variables_clustering_0.35PC.mat'])
    group_model=group{2}; 
    [group_all_all_circle{tk}]=cluster_stability_illustrate_func(time_series,session,neuronIndividuals_new,win_leng,group_model);
end

load('C:\Users\exx\Desktop\HDAC paper fig and method\SFN2019\SFN2019 fig\SCI clustNum plot\tri cir sqr\tri_cir_sqr_dataset.mat');
for tk=1:length(foldername)
    group_all=group_all_all{tk};       
    adjrand_idx_val=[];
    rand_idx_val=[];    
    adjrand_idx_shuffle_val=[];
    rand_idx_shuffle_val=[];    
    for i=1:length(group_all)-1
        [adjrand_idx_val(i),rand_idx_val(i)]=RandIndex(group_all{i},group_all{i+1}); % this calculation is correct, proved by sklearn ari calculation
    end
    [adjrand_idx_val(i+1),rand_idx_val(i+1)]=RandIndex(group_all{3},group_all{1});   
    group_all_mice{tk}=group_all;
    adjrand_idx_val_mice{tk}=adjrand_idx_val;
    rand_idx_val_mice{tk}=rand_idx_val;
end

ori_adjrand=cell2mat(adjrand_idx_val_mice');
ori_rand=cell2mat(rand_idx_val_mice');

% square
for tk=1:length(foldername)
    load([foldername{tk},'\','neuronIndividuals_new.mat']);
    win_leng=5400;%6min
    session=3;
    time_series=1:win_leng/2:size(neuronIndividuals_new{session}.C,2)-win_leng;
    load([foldername{tk},'\','square1_placement_cell_cluster','\','variables_clustering_0.35PC.mat'])
    group_model=group{3}; 
    [group_all_all_square{tk}]=cluster_stability_illustrate_func(time_series,session,neuronIndividuals_new,win_leng,group_model);
end

for tk=1:length(foldername)
    group_all=group_all_all_square{tk};       
    adjrand_idx_val=[];
    rand_idx_val=[];    
    adjrand_idx_shuffle_val=[];
    rand_idx_shuffle_val=[];    
    for i=1:length(group_all)-1
        [adjrand_idx_val(i),rand_idx_val(i)]=RandIndex(group_all{i},group_all{i+1}); % this calculation is correct, proved by sklearn ari calculation
    end
    [adjrand_idx_val(i+1),rand_idx_val(i+1)]=RandIndex(group_all{3},group_all{1});   
    group_all_mice{tk}=group_all;
    adjrand_idx_val_mice{tk}=adjrand_idx_val;
    rand_idx_val_mice{tk}=rand_idx_val;
end

ori_adjrand=cell2mat(adjrand_idx_val_mice');
ori_rand=cell2mat(rand_idx_val_mice');

%% panel N: cluster size at optimal distance, per mice
load('D:\Xu_clusterting_paper_prep11_2020\final_code\final_cluster_data\cluster_region_sz\fig1_3_multiGeo_cluster_Region_sz.mat')
load('D:\Xu_clusterting_paper_prep11_2020\final_code\final_cluster_data\cluster_optimal_num\Fig1_3_multiGeo_clust_original.mat')

circle_sz=[];
circle_shuf_sz=[];
for i=1:6
    circle_sz=[circle_sz;avg_region_multiGeo_k{i,2}{max(group_ori_multiGeo{i,2})}];
    circle_shuf_sz=[circle_shuf_sz;nanmean(avg_region_shuf_multiGeo_k{i,2}{max(group_ori_multiGeo{i,2})})];
end

%% panel O
circle_sz=[];
circle_shuf_sz=[];
for i=1:6
    circle_sz=[circle_sz;avg_region_multiGeo_k{i,3}{max(group_ori_multiGeo{i,3})}];
    circle_shuf_sz=[circle_shuf_sz;nanmean(avg_region_shuf_multiGeo_k{i,3}{max(group_ori_multiGeo{i,3})})];
end
