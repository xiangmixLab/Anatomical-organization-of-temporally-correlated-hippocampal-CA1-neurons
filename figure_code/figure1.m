run('D:\Xu_clusterting_paper_prep11_2020\final code\data_prepare\data_info.m')
all_colors=distinguishable_colors(10);
load(['D:\Xu_clusterting_paper_prep11_2020\final data\final_cluster_data\Fig1_3_multiGeo_clust_original.mat'])

%% A. behav
load([foldername_multiGeo{6},'\','circle1_behav.mat']);
plot(behav.position(:,1),behav.position(:,2),'color',[0.5 0.5 0.5],'lineWidth',2);

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
load([foldername_multiGeo{6},'\','circle1_behav.mat']);

figure;
behav.position=behav.position*280/299;
for i=1:length(segment)
    plot(behav.position(segment{i},1),behav.position(segment{i},2),'color',all_colors(i,:),'lineWidth',2);
    traj_leng(i)=trajectory_length(behav.position(segment{i},:));
    hold on;
end
cen=mean(behav.position,1)-[12 1]; % mean do not directly represent the real center, eyeball adjustment
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

 
%% G: footprint
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
group_period={};
A_color={};
A_color_region={};

for tk=1:length(foldername_multiGeo)
    load([foldername_multiGeo{tk},'\','neuronIndividuals_new.mat']);
    win_leng=5400;%6min
    session=2;
    time_series=1:win_leng/2:size(neuronIndividuals_new{session}.C,2)-win_leng;

    [group_period(tk,:),A_color(tk,:),A_color_region(tk,:)]=cluster_stability_illustrate_func(time_series,session,neuronIndividuals_new,win_leng,group_ori_multiGeo{tk,2});
    
end

for tk=1:size(group_period,2)
    subplot(3,2,2*tk-1)
    imagesc(A_color{6,tk});
    subplot(3,2,2*tk)
    imagesc(A_color_region{6,tk});
end

%% I: region size, all 6 mice 
avg_max_region={};
avg_max_reg_shuf_all={};
for tk=1:6
    for k=2:10
        load([foldername_multiGeo{tk},'\neuronIndividuals_new.mat']);
        clust_num=[2:10];
        shuffles_num=100;
        cond_sign=2;

        [~,~,avg_max_region{tk}(k),avg_max_reg_shuf_all{tk}(k)]=DBSCAN_region_quantify_022422(group_ori_multiGeo{tk,2},neuronIndividuals_new,[]);
    end
end

for tk=1:6
   avg_max_region_mat(tk,:)=avg_max_region{tk}(2,:);
   avg_max_reg_shuf_all_mat(tk,:)=avg_max_reg_shuf_all{tk}(2,:);
end

mean_avg_max_region_mat=mean(avg_max_region_mat,1);
se_avg_max_region_mat=std(avg_max_region_mat,[],1)/(size(avg_max_region_mat,1)^0.5);
min_avg_max_region_mat=mean(avg_max_region_mat,1)-se_avg_max_region_mat;
max_avg_max_region_mat=mean(avg_max_region_mat,1)+se_avg_max_region_mat;

mean_avg_max_reg_shuf_all_mat=mean(avg_max_reg_shuf_all_mat,1);
se_avg_max_reg_shuf_all_mat=std(avg_max_reg_shuf_all_mat,[],1)/(size(avg_max_reg_shuf_all_mat,1)^0.5);
min_avg_max_reg_shuf_all_mat=mean(avg_max_reg_shuf_all_mat,1)-se_avg_max_reg_shuf_all_mat;
max_avg_max_reg_shuf_all_mat=mean(avg_max_reg_shuf_all_mat,1)+se_avg_max_reg_shuf_all_mat;

x=[[0:length(mean_avg_max_region_mat)-1]',[0:length(mean_avg_max_reg_shuf_all_mat)-1]'];
y=[mean_avg_max_region_mat',mean_avg_max_reg_shuf_all_mat'];
minMaxRange={[min_avg_max_region_mat',max_avg_max_region_mat'],[min_avg_max_reg_shuf_all_mat',max_avg_max_reg_shuf_all_mat']};
% plotwithMaxMixRange(x,y,minMaxRange,{[0 0 0.9],[0 0 0.55],[0 0 0.20]},{'-o','-s','-^'});

errorbar(0:length(mean_avg_max_region_mat)-1,mean_avg_max_region_mat,se_avg_max_region_mat,'-o','color',[0    0.4470    0.7410]);
hold on;
errorbar(0:length(mean_avg_max_reg_shuf_all_mat)-1,mean_avg_max_reg_shuf_all_mat,se_avg_max_reg_shuf_all_mat,'-s','color',[0.8500    0.3250    0.0980]);

xlim([-1 9])

%% J: distance-temporal correlation plot
cluster_file_name={};
for i=1:6
    fname=foldername{mice_all(tk)};
    cluster_file_name{tk}=[fname,'\','circle1_placement_cell_cluster','\','variables_clustering_0.35PC.mat'];
    session=2;
end
[all_dis_corr_circle_all_intra,all_dis_corr_circle_all_inter,all_dis_corr_circle_all_all,midx,f,gof]=pairwise_dis_tempCorr_091421(foldername,cluster_file_name,session);

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
%% panel h cross time stability
group_all_mice={};
adjoverlap_val_mice={};
overlap_val_mice={};
adjoverlap_val_shuffle_mice={};
overlap_val_shuffle_mice={};
adjoverlap_val_shuffle_95_mice={};
overlap_val_shuffle_95_mice={};
overlap_val_shuffle_nc_mice={};

for tk=1:length(foldername)
    load([foldername{tk},'\','neuronIndividuals_new.mat']);
    win_leng=5400;%6min
    session=2;
    time_series=1:win_leng/2:size(neuronIndividuals_new{session}.C,2)-win_leng;

    load([foldername{tk},'\','circle1_placement_cell_cluster','\','variables_clustering_0.35PC.mat'])
    group_curr=group{2}; 

    [group_all_all{tk}]=cluster_stability_illustrate_func(time_series,session,neuronIndividuals_new,win_leng,group_curr);
end

for tk=1:length(foldername)
    group_all=group_all_all{tk};
    figure;
    load([foldername{tk},'\','neuronIndividuals_new.mat']);

    for i=1:length(group_all)
        subplot(2,9,i);
        [A_color,A_color_region,avg_max_region,avg_max_reg_shuf_all,avg_max_reg_shuf_all_95]=DBSCAN_region_quantify_func_single_group(group_all{i},2,neuronIndividuals_new);
        imagesc(A_color);
        subplot(2,9,i+9);
        imagesc(A_color_region);
    end
    
    group_all_mice{tk}=group_all;
end

for tk=1:length(foldername)
    group_all=group_all_mice{tk};
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


% ori_adjrand=cell2mat(adjoverlap_val_mice');
ori_rand=cell2mat(overlap_val_mice');
% ori_rand_nc=cell2mat(overlap_val_shuffle_nc_mice');

% ori_adjrand_shuffle_qt=cell2mat(adjoverlap_val_shuffle_qt_mice)';
ori_rand_shuffle_qt=cell2mat(overlap_val_shuffle_qt_mice)';

rand_diff_from_shuffle=ori_rand-ori_rand_shuffle_qt;
% shuffle_adjrand_mean=[];
% shuffle_adjrand_95=[];
% for i=1:length(foldername)   
%     mean_adjrand_val_shuffle_mice=mean(adjoverlap_val_shuffle_mice{i},1);
%     adjrand_val_shuffle_95_mice=quantile(adjoverlap_val_shuffle_mice{i},.95,1);
%     shuffle_adjrand_mean(i,:)=mean_adjrand_val_shuffle_mice;
%     shuffle_adjrand_95(i,:)=adjrand_val_shuffle_95_mice;
% end   
% max(max(shuffle_adjrand_95))    

% overtime significance
ctt=1;
p_value_list=[];
for i=1:size(ori_adjrand,2)-1
    for j=i+1:size(ori_adjrand,2)
        p_value_list(ctt)=ranksum(ori_adjrand(:,i),ori_adjrand(:,j));
        ctt=ctt+1;
    end
end
p_value_mat=squareform(p_value_list);

%% Panel c: autocorrelation of cluster temporal trace
for j = 1:length(unique(group_curr))
    rangee=find(group_curr == j);
    nC=neuronIndividuals_new{2}.C(rangee,1:end)';
    nC_sum=sum(nC,1);
    [acf{j},lag{j}] = xcorr(nC_sum','coeff');
end

for j = 1:length(unique(group_curr))
    subplot(length(unique(group_curr)),1,j);
    plot(lag{j},acf{j},'-','color',colorClusters2(j,:));hold on;
    xlim([-75 75])
    ctt=ctt+1;
end

%% panel n: similarity over time
load('D:\Remapping_square_circle_triangle_061119_061319\M3425F\group_all')

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
    cond_sign=2;

    load([foldername{tk},'\circle1_placement_cell_cluster','\','variables_clustering_0.35PC.mat'])
    group_curr=group{2};
    clust_idx=length(unique(group_curr));
    gp_rec{cond_sign}{clust_idx-1}=group_curr;

    [~,~,avg_max_region{tk},avg_max_reg_shuf_all{tk},avg_max_reg_shuf_all_95{tk}]=DBSCAN_region_quantify_func_simplify_no_plot(gp_rec,clust_num,shuffles_num,clust_idx-1,cond_sign,neuronIndividuals_new);


    avg_max_region_mat=avg_max_region{tk}(2,:);
    avg_max_reg_shuf_all_mat=avg_max_reg_shuf_all{tk}(2,:);
    avg_max_reg_shuf_all_95_mat=avg_max_reg_shuf_all_95{tk}(2,:);

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
end



 %% new D : average dis in 4.583 sec, circle1
t_leng1=30*m_hlhh_circle; % for all mice
t_leng2=30*m_hlhh_square; % for all mice

distance_sec_cir_mean=[];
distance_sec_sqr_mean=[];
cir_radius=zeros(1,length(foldername));
sqr_radius=zeros(1,length(foldername));
for i=1:length(foldername)
    load([foldername{i},'\','circle1_Behav.mat']);
    dis_all_sec=[];
    for j=1:size(behav.position,1)-t_leng1(i)
        position_sec=behav.position(j:j+t_leng1(i)-1,:);
        distance_sec=0;
        for k=1:size(position_sec,1)-1
            distance_sec=distance_sec+nansum((position_sec(k+1,:)-position_sec(k,:)).^2).^0.5;
        end
        dis_sec_all(j)=distance_sec;
    end
    distance_sec_cir_mean(i)=nanmean(dis_sec_all);
    
    cir_radius(i)=mean([max(behav.position(:,1))-min(behav.position(:,1)),max(behav.position(:,2))-min(behav.position(:,2))]);

    
    load([foldername{i},'\','square1_Behav.mat']);
    dis_sec_all=[];
    for j=1:size(behav.position,1)-t_leng2(i)
        position_sec=behav.position(j:j+t_leng2(i)-1,:);
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

%% new f: magnification of 50sec segment
