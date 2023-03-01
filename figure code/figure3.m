%% Figure 3: Anatomical cluster-specific calcium event activities vary across different environments.

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% OF secion %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
foldername={
    'D:\Xu_clusterting_paper_prep11_2020\final data\final_neuron_behav_data\Fig_1_3_tri_cir_sqr\M3411\';
    'D:\Xu_clusterting_paper_prep11_2020\final data\final_neuron_behav_data\Fig_1_3_tri_cir_sqr\M3412\';
    'D:\Xu_clusterting_paper_prep11_2020\final data\final_neuron_behav_data\Fig_1_3_tri_cir_sqr\M3421F\';
    'D:\Xu_clusterting_paper_prep11_2020\final data\final_neuron_behav_data\Fig_1_3_tri_cir_sqr\M3422F\';
    'D:\Xu_clusterting_paper_prep11_2020\final data\final_neuron_behav_data\Fig_1_3_tri_cir_sqr\M3424F\';
    'D:\Xu_clusterting_paper_prep11_2020\final data\final_neuron_behav_data\Fig_1_3_tri_cir_sqr\M3425F\'
    }
behavName={
    'triangle1_Behav.mat';
    'circle1_Behav.mat';
    'square1_Behav.mat';    
    'circle2_Behav.mat';
    'square2_Behav.mat';
    'triangle2_Behav.mat';
    }

cluster_filename='D:\Xu_clusterting_paper_prep11_2020\final data\final_cluster_data\Fig1_3_multiGeo_clust_original.mat'
mice=2; % chosen example mice for behav and footprint illustration

%% load all behav data
all_behav={};
for i=1:length(foldername)
    for j=1:length(behavName)
        load([foldername{1},'\',behavName{j}]);
        all_behav{i,j}=behav;
    end
end

%% A. OF traj
figure;
for i=1:length(behavName)
    subplot(2,3,i);
    plot(all_behav{mice,i}.position(:,1),all_behav{mice,i}.position(:,2),'color',[0.5 0.5 0.5],'lineWidth',2);
    set(gcf,'renderer','painters');
end

%% B: tri-cir-sqr footprint

load([foldername{mice},'\','neuronIndividuals_new.mat'])
load(cluster_filename)

group=group_ori_multiGeo_sameNum_eachMice(mice,:);
for j=1:length(neuronIndividuals_new)  
    [~,group_cond_j_reindexed]=cluster_reindex_112522(group{1},group{j}); % try to reIndex the clusters to make the color assignment similar to the first trial
    [A_color{1,j},A_color_region{1,j}]=DBSCAN_region_quantify_022422(group_cond_j_reindexed,neuronIndividuals_new,[]);
end

% plot footprint
figure;
for j=1:length(neuronIndividuals_new)  
    subplot(2,3,j);
    imagesc(imrotate(A_color{1,j},-135));
end

% plot region
figure;
for j=1:length(neuronIndividuals_new)  
    subplot(2,3,j);
    imagesc(imrotate(A_color_region{1,j},-135));
end


%% C: pearson correlation for first/second half trials, same geo
load('D:\final_HDAC_AD_automatic_processing\dataloc\Remapping_square_circle_triangle_061119_061319_merged.mat')

rateMap_corr_all_halves_of=[];
rateMap_corr_all_halves_of_t=[];
rateMap_corr_all_shuf_halves_of=[];
tic;
load(['Y:\Lujia\sfn2019\SFN2019 fig and text\Round15\fig3 panels\of_temporal_half_data.mat']);
for tk=1:6
   rateMap_corr1=rateMap_correlation(behavpos_half_ratemap_of{tk,1}{1},behavpos_half_ratemap_of{tk,1}{2},behavpos_half_ct_of{tk,1}{1},behavpos_half_ct_of{tk,1}{2},1);
   rateMap_corr2=rateMap_correlation(behavpos_half_ratemap_of{tk,2}{1},behavpos_half_ratemap_of{tk,2}{2},behavpos_half_ct_of{tk,2}{1},behavpos_half_ct_of{tk,2}{2},1);
   rateMap_corr3=rateMap_correlation(behavpos_half_ratemap_of{tk,3}{1},behavpos_half_ratemap_of{tk,3}{2},behavpos_half_ct_of{tk,3}{1},behavpos_half_ct_of{tk,3}{2},1);
   rateMap_corr4=rateMap_correlation(behavpos_half_ratemap_of{tk,4}{1},behavpos_half_ratemap_of{tk,4}{2},behavpos_half_ct_of{tk,4}{1},behavpos_half_ct_of{tk,4}{2},1);
   rateMap_corr5=rateMap_correlation(behavpos_half_ratemap_of{tk,5}{1},behavpos_half_ratemap_of{tk,5}{2},behavpos_half_ct_of{tk,5}{1},behavpos_half_ct_of{tk,5}{2},1);
   rateMap_corr6=rateMap_correlation(behavpos_half_ratemap_of{tk,6}{1},behavpos_half_ratemap_of{tk,6}{2},behavpos_half_ct_of{tk,6}{1},behavpos_half_ct_of{tk,6}{2},1);

   rateMap_corr_all_halves_of_t(tk,:)=[[nanmean(rateMap_corr1),nanmean(rateMap_corr2),nanmean(rateMap_corr3),nanmean(rateMap_corr4),nanmean(rateMap_corr5),nanmean(rateMap_corr6)]];
   rateMap_corr_all_halves_of(tk,1)=mean([nanmean(rateMap_corr1),nanmean(rateMap_corr2),nanmean(rateMap_corr3),nanmean(rateMap_corr4),nanmean(rateMap_corr5),nanmean(rateMap_corr6)]);

   parfor i=1:1000
        idx_rand_fr1=randperm(length(behavpos_half_ratemap_of{tk,1}{1}));
        idx_rand_fr2=randperm(length(behavpos_half_ratemap_of{tk,1}{2}));
        idx_rand_fr3=randperm(length(behavpos_half_ratemap_of{tk,2}{1}));
        idx_rand_fr4=randperm(length(behavpos_half_ratemap_of{tk,2}{2}));
        idx_rand_fr5=randperm(length(behavpos_half_ratemap_of{tk,3}{1}));
        idx_rand_fr6=randperm(length(behavpos_half_ratemap_of{tk,3}{2}));
        idx_rand_fr7=randperm(length(behavpos_half_ratemap_of{tk,4}{1}));
        idx_rand_fr8=randperm(length(behavpos_half_ratemap_of{tk,4}{2}));
        idx_rand_fr9=randperm(length(behavpos_half_ratemap_of{tk,5}{1}));
        idx_rand_fr10=randperm(length(behavpos_half_ratemap_of{tk,5}{2}));
        idx_rand_fr11=randperm(length(behavpos_half_ratemap_of{tk,6}{1}));
        idx_rand_fr12=randperm(length(behavpos_half_ratemap_of{tk,6}{2}));
        
        fr1_shuf=behavpos_half_ratemap_of{tk,1}{1}(idx_rand_fr1);
        ct1_shuf=behavpos_half_ct_of{tk,1}{1};
        fr2_shuf=behavpos_half_ratemap_of{tk,1}{2}(idx_rand_fr2);
        ct2_shuf=behavpos_half_ct_of{tk,1}{2};
        fr3_shuf=behavpos_half_ratemap_of{tk,2}{1}(idx_rand_fr3);
        ct3_shuf=behavpos_half_ct_of{tk,2}{1}; 
        fr4_shuf=behavpos_half_ratemap_of{tk,2}{2}(idx_rand_fr4);
        ct4_shuf=behavpos_half_ct_of{tk,2}{2};
        fr5_shuf=behavpos_half_ratemap_of{tk,3}{1}(idx_rand_fr5);
        ct5_shuf=behavpos_half_ct_of{tk,3}{1};
        fr6_shuf=behavpos_half_ratemap_of{tk,3}{2}(idx_rand_fr6);
        ct6_shuf=behavpos_half_ct_of{tk,3}{2};
        fr7_shuf=behavpos_half_ratemap_of{tk,4}{1}(idx_rand_fr7);
        ct7_shuf=behavpos_half_ct_of{tk,4}{1};
        fr8_shuf=behavpos_half_ratemap_of{tk,4}{2}(idx_rand_fr8);
        ct8_shuf=behavpos_half_ct_of{tk,4}{2};
        fr9_shuf=behavpos_half_ratemap_of{tk,5}{1}(idx_rand_fr9);
        ct9_shuf=behavpos_half_ct_of{tk,5}{1}; 
        fr10_shuf=behavpos_half_ratemap_of{tk,5}{2}(idx_rand_fr10);
        ct10_shuf=behavpos_half_ct_of{tk,5}{2};
        fr11_shuf=behavpos_half_ratemap_of{tk,6}{1}(idx_rand_fr11);
        ct11_shuf=behavpos_half_ct_of{tk,6}{1};
        fr12_shuf=behavpos_half_ratemap_of{tk,6}{2}(idx_rand_fr12);
        ct12_shuf=behavpos_half_ct_of{tk,6}{2};
        
       rateMap_corr1_shuf=rateMap_correlation(fr1_shuf,fr2_shuf,ct1_shuf,ct2_shuf,1);
       rateMap_corr2_shuf=rateMap_correlation(fr3_shuf,fr4_shuf,ct3_shuf,ct4_shuf,1);
       rateMap_corr3_shuf=rateMap_correlation(fr5_shuf,fr6_shuf,ct5_shuf,ct6_shuf,1);
       rateMap_corr4_shuf=rateMap_correlation(fr7_shuf,fr8_shuf,ct7_shuf,ct8_shuf,1);
       rateMap_corr5_shuf=rateMap_correlation(fr9_shuf,fr10_shuf,ct9_shuf,ct10_shuf,1);
       rateMap_corr6_shuf=rateMap_correlation(fr11_shuf,fr12_shuf,ct11_shuf,ct12_shuf,1);
       rateMap_corr_all_shuf_halves_of_t(tk,:,i)=[nanmean(rateMap_corr1_shuf),nanmean(rateMap_corr2_shuf),nanmean(rateMap_corr3_shuf),nanmean(rateMap_corr4_shuf),nanmean(rateMap_corr5_shuf),nanmean(rateMap_corr6_shuf)];       
       rateMap_corr_all_shuf_halves_of(tk,1,i)=mean([nanmean(rateMap_corr1_shuf),nanmean(rateMap_corr2_shuf),nanmean(rateMap_corr3_shuf),nanmean(rateMap_corr4_shuf),nanmean(rateMap_corr5_shuf),nanmean(rateMap_corr6_shuf)]);
    end
toc;
end

rateMap_corr_all_shuf_halves_of1=quantile(rateMap_corr_all_shuf_halves_of,.95,3);

shuffle=1000;
rateMap_corr_all_of_t=[];
rateMap_corr_all_of=[];
rateMap_corr_all_of_cs=[];
rateMap_corr_all_shuf_of=[];
rateMap_corr_all_shuf_of_cs=[];
rateMap_corr_all_shuf_of_t=[];
tic;
for tk=1:6
    load([foldername{tk},'\neuronIndividuals_new.mat']);
    
   [fr_all,ct_all]=figure3_loadFiringRate_G(foldername{tk});
    
   rateMap_corr1=rateMap_correlation(fr_all{1},fr_all{6},ct_all{1},ct_all{6},1);
   rateMap_corr2=rateMap_correlation(fr_all{2},fr_all{4},ct_all{2},ct_all{4},1);  
   rateMap_corr3=rateMap_correlation(fr_all{3},fr_all{5},ct_all{3},ct_all{5},1);  
   rateMap_corr4=rateMap_correlation(fr_all{1},fr_all{2},ct_all{1},ct_all{2},1);
   rateMap_corr5=rateMap_correlation(fr_all{2},fr_all{3},ct_all{2},ct_all{3},1);  
   rateMap_corr6=rateMap_correlation(fr_all{1},fr_all{3},ct_all{1},ct_all{3},1);  
   rateMap_corr7=rateMap_correlation(fr_all{6},fr_all{4},ct_all{6},ct_all{4},1);
   rateMap_corr8=rateMap_correlation(fr_all{4},fr_all{5},ct_all{4},ct_all{5},1);  
   rateMap_corr9=rateMap_correlation(fr_all{6},fr_all{5},ct_all{6},ct_all{5},1);
   rateMap_corr10=rateMap_correlation(fr_all{1},fr_all{4},ct_all{1},ct_all{4},1);
   rateMap_corr11=rateMap_correlation(fr_all{1},fr_all{5},ct_all{1},ct_all{5},1);  
   rateMap_corr12=rateMap_correlation(fr_all{2},fr_all{5},ct_all{2},ct_all{5},1);
   rateMap_corr13=rateMap_correlation(fr_all{2},fr_all{6},ct_all{2},ct_all{6},1);
   rateMap_corr14=rateMap_correlation(fr_all{3},fr_all{4},ct_all{3},ct_all{4},1);
   rateMap_corr15=rateMap_correlation(fr_all{3},fr_all{6},ct_all{3},ct_all{6},1);
   rateMap_corr_all_of_t(tk,:)=[nanmean(rateMap_corr1),nanmean(rateMap_corr2),nanmean(rateMap_corr3),nanmean(rateMap_corr4),nanmean(rateMap_corr5),nanmean(rateMap_corr6),nanmean(rateMap_corr7),nanmean(rateMap_corr8),nanmean(rateMap_corr9),nanmean(rateMap_corr10),nanmean(rateMap_corr11),nanmean(rateMap_corr12),nanmean(rateMap_corr13),nanmean(rateMap_corr14),nanmean(rateMap_corr15)];
   rateMap_corr_all_of(tk,:)=[mean([nanmean(rateMap_corr1),nanmean(rateMap_corr2),nanmean(rateMap_corr3)],2),mean([nanmean(rateMap_corr4),nanmean(rateMap_corr5),nanmean(rateMap_corr6),nanmean(rateMap_corr7),nanmean(rateMap_corr8),nanmean(rateMap_corr9)],2),mean([nanmean(rateMap_corr10),nanmean(rateMap_corr11),nanmean(rateMap_corr12),nanmean(rateMap_corr13),nanmean(rateMap_corr14),nanmean(rateMap_corr15)],2)];
   rateMap_corr_all_of_cs(tk,:)=[mean([nanmean(rateMap_corr4),nanmean(rateMap_corr5),nanmean(rateMap_corr6),nanmean(rateMap_corr7),nanmean(rateMap_corr8),nanmean(rateMap_corr9)],2),mean([nanmean(rateMap_corr10),nanmean(rateMap_corr11),nanmean(rateMap_corr12),nanmean(rateMap_corr13),nanmean(rateMap_corr14),nanmean(rateMap_corr15)],2)];

   %      f1=fit(rateMap_corr',temp_corr','power2');
    for j=1:shuffle
        idx_fr1=randperm(length(fr_all{1}));
        idx_fr2=randperm(length(fr_all{2}));
        idx_fr3=randperm(length(fr_all{3}));
        idx_fr4=randperm(length(fr_all{4}));
        idx_fr5=randperm(length(fr_all{5}));
        idx_fr6=randperm(length(fr_all{6}));
        
        fr_all1=fr_all{1}(idx_fr1);
        ct_all1=ct_all{1};
        fr_all2=fr_all{2}(idx_fr2);
        ct_all2=ct_all{2};        
        fr_all3=fr_all{3}(idx_fr3);
        ct_all3=ct_all{3};       
        fr_all4=fr_all{4}(idx_fr4);
        ct_all4=ct_all{4}; 
        fr_all5=fr_all{5}(idx_fr5);
        ct_all5=ct_all{5};       
        fr_all6=fr_all{6}(idx_fr6);
        ct_all6=ct_all{6};
        
        rateMap_corr1=rateMap_correlation(fr_all1,fr_all6,ct_all1,ct_all6,1);
        rateMap_corr2=rateMap_correlation(fr_all2,fr_all4,ct_all2,ct_all4,1);  
        rateMap_corr3=rateMap_correlation(fr_all3,fr_all5,ct_all3,ct_all5,1);  
        rateMap_corr4=rateMap_correlation(fr_all1,fr_all2,ct_all1,ct_all2,1);
        rateMap_corr5=rateMap_correlation(fr_all2,fr_all3,ct_all2,ct_all3,1);  
        rateMap_corr6=rateMap_correlation(fr_all1,fr_all3,ct_all1,ct_all3,1);  
        rateMap_corr7=rateMap_correlation(fr_all6,fr_all4,ct_all6,ct_all4,1);
        rateMap_corr8=rateMap_correlation(fr_all4,fr_all5,ct_all4,ct_all5,1);  
        rateMap_corr9=rateMap_correlation(fr_all6,fr_all5,ct_all6,ct_all5,1);
        rateMap_corr10=rateMap_correlation(fr_all1,fr_all4,ct_all1,ct_all4,1);
        rateMap_corr11=rateMap_correlation(fr_all1,fr_all5,ct_all1,ct_all5,1);  
        rateMap_corr12=rateMap_correlation(fr_all2,fr_all5,ct_all2,ct_all5,1);
        rateMap_corr13=rateMap_correlation(fr_all2,fr_all6,ct_all2,ct_all6,1);
        rateMap_corr14=rateMap_correlation(fr_all3,fr_all4,ct_all3,ct_all4,1);
        rateMap_corr15=rateMap_correlation(fr_all3,fr_all6,ct_all3,ct_all6,1);
        rateMap_corr_all_shuf_of_t(tk,:,j)=[nanmean(rateMap_corr1),nanmean(rateMap_corr2),nanmean(rateMap_corr3),nanmean(rateMap_corr4),nanmean(rateMap_corr5),nanmean(rateMap_corr6),nanmean(rateMap_corr7),nanmean(rateMap_corr8),nanmean(rateMap_corr9),nanmean(rateMap_corr10),nanmean(rateMap_corr11),nanmean(rateMap_corr12),nanmean(rateMap_corr13),nanmean(rateMap_corr14),nanmean(rateMap_corr15)];
        rateMap_corr_all_shuf_of(tk,:,j)=[mean([nanmean(rateMap_corr1),nanmean(rateMap_corr2),nanmean(rateMap_corr3)],2),mean([nanmean(rateMap_corr4),nanmean(rateMap_corr5),nanmean(rateMap_corr6),nanmean(rateMap_corr7),nanmean(rateMap_corr8),nanmean(rateMap_corr9)],2),mean([nanmean(rateMap_corr10),nanmean(rateMap_corr11),nanmean(rateMap_corr12),nanmean(rateMap_corr13),nanmean(rateMap_corr14),nanmean(rateMap_corr15)],2)];
        rateMap_corr_all_shuf_of_cs(tk,:,j)=[mean([nanmean(rateMap_corr4),nanmean(rateMap_corr5),nanmean(rateMap_corr6),nanmean(rateMap_corr7),nanmean(rateMap_corr8),nanmean(rateMap_corr9)],2),mean([nanmean(rateMap_corr10),nanmean(rateMap_corr11),nanmean(rateMap_corr12),nanmean(rateMap_corr13),nanmean(rateMap_corr14),nanmean(rateMap_corr15)],2)];
    end
toc;
end
rateMap_corr_all_of=[rateMap_corr_all_halves_of,rateMap_corr_all_of];
rateMap_corr_all_shuf_of1=quantile(rateMap_corr_all_shuf_of,.95,3);
rateMap_corr_all_shuf_of1_t=quantile(rateMap_corr_all_shuf_of_t,.95,3);
rateMap_corr_all_shuf_of_cs=quantile(rateMap_corr_all_shuf_of_cs,.95,3);

%% c: cluster overlap OF

shuffle=1000;
tri1_tri2_rand=[];
cir1_cir2_rand=[];
sqr1_sqr2_rand=[];
tri1_cir1_rand=[];
cir1_sqr1_rand=[];
tir1_sqr1_rand=[];
tri2_cir2_rand=[];
cir2_sqr2_rand=[];
tri2_sqr2_rand=[];
tri1_cir2_rand=[];
tri1_sqr2_rand=[];
cir1_tri2_rand=[];
cir1_sqr2_rand=[];
sqr1_tri2_rand=[];
sqr1_cir2_rand=[];

tri1_tri2_adjrand=[];
cir1_cir2_adjrand=[];
sqr1_sqr2_adjrand=[];
tri1_cir1_adjrand=[];
cir1_sqr1_adjrand=[];
tir1_sqr1_adjrand=[];
tri2_cir2_adjrand=[];
cir2_sqr2_adjrand=[];
tri2_sqr2_adjrand=[];
tri1_cir2_adjrand=[];
tri1_sqr2_adjrand=[];
cir1_tri2_adjrand=[];
cir1_sqr2_adjrand=[];
sqr1_tri2_adjrand=[];
sqr1_cir2_adjrand=[];

tri1_tri2_shuffle_rand=[];
cir1_cir2_shuffle_rand=[];
sqr1_sqr2_shuffle_rand=[];
tri1_cir1_shuffle_rand=[];
cir1_sqr1_shuffle_rand=[];
tir1_sqr1_shuffle_rand=[];
tri2_cir2_shuffle_rand=[];
cir2_sqr2_shuffle_rand=[];
tri2_sqr2_shuffle_rand=[];
tri1_cir2_shuffle_rand=[];
tri1_sqr2_shuffle_rand=[];
cir1_tri2_shuffle_rand=[];
cir1_sqr2_shuffle_rand=[];
sqr1_tri2_shuffle_rand=[];
sqr1_cir2_shuffle_rand=[];

clust_idx=4;

mice_all=[1:6]

for i=1:length(mice_all)
    mice=mice_all(i);
    group=group_record(:,mice);
    g1=group{1}{clust_idx};
    g2=group{2}{clust_idx};
    g3=group{3}{clust_idx};
    g4=group{4}{clust_idx};
    g5=group{5}{clust_idx};
    g6=group{6}{clust_idx};

    [tri1_tri2_rand(i,1)]=new_cluster_overlap(g1,g6);
    [cir1_cir2_rand(i,1)]=new_cluster_overlap(g2,g4);
    [sqr1_sqr2_rand(i,1)]=new_cluster_overlap(g3,g5);
    [tri1_cir1_rand(i,1)]=new_cluster_overlap(g1,g2);
    [cir1_sqr1_rand(i,1)]=new_cluster_overlap(g2,g3);
    [tir1_sqr1_rand(i,1)]=new_cluster_overlap(g1,g3);
    [tri2_cir2_rand(i,1)]=new_cluster_overlap(g6,g4);
    [cir2_sqr2_rand(i,1)]=new_cluster_overlap(g4,g5);
    [tri2_sqr2_rand(i,1)]=new_cluster_overlap(g6,g5);
    [tri1_cir2_rand(i,1)]=new_cluster_overlap(g1,g4);
    [tri1_sqr2_rand(i,1)]=new_cluster_overlap(g1,g5);
    [cir1_tri2_rand(i,1)]=new_cluster_overlap(g2,g5);
    [cir1_sqr2_rand(i,1)]=new_cluster_overlap(g2,g6);
    [sqr1_tri2_rand(i,1)]=new_cluster_overlap(g3,g4);
    [sqr1_cir2_rand(i,1)]=new_cluster_overlap(g3,g6);

    
    for j=1:shuffle
        g1=group{1}{clust_idx}(randperm(length(group{1}{clust_idx})));
        g2=group{2}{clust_idx}(randperm(length(group{2}{clust_idx})));
        g3=group{3}{clust_idx}(randperm(length(group{3}{clust_idx})));
        g4=group{4}{clust_idx}(randperm(length(group{4}{clust_idx})));
        g5=group{5}{clust_idx}(randperm(length(group{5}{clust_idx})));
        g6=group{6}{clust_idx}(randperm(length(group{6}{clust_idx})));       

        [tri1_tri2_shuffle_rand(i,j)]=new_cluster_overlap(g1,g6);
        [cir1_cir2_shuffle_rand(i,j)]=new_cluster_overlap(g2,g4);
        [sqr1_sqr2_shuffle_rand(i,j)]=new_cluster_overlap(g3,g5);
        [tri1_cir1_shuffle_rand(i,j)]=new_cluster_overlap(g1,g2);
        [cir1_sqr1_shuffle_rand(i,j)]=new_cluster_overlap(g2,g3);
        [tir1_sqr1_shuffle_rand(i,j)]=new_cluster_overlap(g1,g3);
        [tri2_cir2_shuffle_rand(i,j)]=new_cluster_overlap(g6,g4);
        [cir2_sqr2_shuffle_rand(i,j)]=new_cluster_overlap(g4,g5);
        [tri2_sqr2_shuffle_rand(i,j)]=new_cluster_overlap(g6,g5);
        [tri1_cir2_shuffle_rand(i,j)]=new_cluster_overlap(g1,g4);
        [tri1_sqr2_shuffle_rand(i,j)]=new_cluster_overlap(g1,g5);
        [cir1_tri2_shuffle_rand(i,j)]=new_cluster_overlap(g2,g5);
        [cir1_sqr2_shuffle_rand(i,j)]=new_cluster_overlap(g2,g6);
        [sqr1_tri2_shuffle_rand(i,j)]=new_cluster_overlap(g3,g4);
        [sqr1_cir2_shuffle_rand(i,j)]=new_cluster_overlap(g3,g6);
    end
end
rand_all=[mean([tri1_tri2_rand,cir1_cir2_rand,sqr1_sqr2_rand],2),mean([tri1_cir1_rand,cir1_sqr1_rand,tir1_sqr1_rand,tri2_cir2_rand,cir2_sqr2_rand,tri2_sqr2_rand],2),mean([tri1_cir2_rand,tri1_sqr2_rand,cir1_sqr2_rand,cir1_tri2_rand,sqr1_cir2_rand,sqr1_tri2_rand],2)];
% rand_adjall=[mean([tri1_tri2_adjrand,cir1_cir2_adjrand,sqr1_sqr2_adjrand],2),mean([tri1_cir1_adjrand,cir1_sqr1_adjrand,tir1_sqr1_adjrand,tri2_cir2_adjrand,cir2_sqr2_adjrand,tri2_sqr2_adjrand],2),mean([tri1_cir2_adjrand,tri1_sqr2_adjrand,cir1_sqr2_adjrand,cir1_tri2_adjrand,sqr1_cir2_adjrand,sqr1_tri2_adjrand],2)];

dat={tri1_tri2_shuffle_rand;cir1_cir2_shuffle_rand;sqr1_sqr2_shuffle_rand;tri1_cir1_shuffle_rand;cir1_sqr1_shuffle_rand;tir1_sqr1_shuffle_rand;tri2_cir2_shuffle_rand;cir2_sqr2_shuffle_rand;tri2_sqr2_shuffle_rand;tri1_cir2_shuffle_rand;tri1_sqr2_shuffle_rand;cir1_tri2_shuffle_rand;cir1_sqr2_shuffle_rand;sqr1_tri2_shuffle_rand;sqr1_cir2_shuffle_rand};
dat_95=[];
for i=1:length(dat)
    dat_95(:,i)=quantile(dat{i},.95,2);
end
dat_95_mean=[mean(dat_95(:,1:3),2),mean(dat_95(:,4:9),2),mean(dat_95(:,10:15),2)];

ranksum_test_summary=[];
for i=1:size(rand_all,2)-1
    for j=i+1:size(rand_all,2)
        ranksum_test_summary(i,j)=ranksum(rand_all(:,i),rand_all(:,j));
    end
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% linear track section %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('D:\final_HDAC_AD_automatic_processing\dataloc\Remapping_linear_track_053019_053119_merged.mat')
load('C:\Users\exx\Desktop\HDAC paper fig and method\SFN2019\SFN2019 fig\SCI clustNum plot\linear track\dataset_water_removed.mat');
foldername=unique(destination);

%% d1: linear track behav track
figure;
load([foldername{1},'\','neuronIndividuals_new.mat'])
subplot(131);
load([foldername{1},'\','Horizontal1_Behav.mat'])
behavpos_new=behav.position;
plot(behavpos_new(:,1),behavpos_new(:,2),'color',[0.5 0.5 0.5],'lineWidth',2);
xlim([0 1000]);
subplot(132);
load([foldername{1},'\','Horizontal2_Behav.mat'])
behavpos_new=behav.position;
plot(behavpos_new(:,1),behavpos_new(:,2),'color',[0.5 0.5 0.5],'lineWidth',2);
xlim([0 1000]);
subplot(133);
load([foldername{1},'\','Vertical_Behav.mat'])
behavpos_new=behav.position;
plot(behavpos_new(:,1),behavpos_new(:,2),'color',[0.5 0.5 0.5],'lineWidth',2);
ylim([0 1000]);


%% d2
mice_all=[4];
for i=1:length(mice_all)
    mice=mice_all(i)
    load([foldername{mice},'\','neuronIndividuals_new.mat'])
    group=group_record(:,mice);
    for j=1:length(neuronIndividuals_new)
        mice=mice_all(i);
        clust_num=[2:10];
        shuffles_num=100;
        clust_idx=4;
        cond_sign=j;
        [~,group1in2]=determineSharedCells_new(group{1}{clust_idx},group{cond_sign}{clust_idx});
        group{cond_sign}{clust_idx}=group1in2;
        DBSCAN_region_quantify_func_simplify(group,clust_num,shuffles_num,clust_idx,cond_sign,neuronIndividuals_new)
        figure(1)
        saveas(gcf,['C:\Users\exx\Desktop\HDAC paper fig and method\SFN2019\SFN2019 fig\Round9\fig3 panels\b',num2str(j),'1.eps'],'epsc');
        figure(3)
        saveas(gcf,['C:\Users\exx\Desktop\HDAC paper fig and method\SFN2019\SFN2019 fig\Round9\fig3 panels\b',num2str(j),'2.eps'],'epsc');
        figure(4)
        saveas(gcf,['C:\Users\exx\Desktop\HDAC paper fig and method\SFN2019\SFN2019 fig\Round9\fig3 panels\b',num2str(j),'3.eps'],'epsc');
        close all
    end
end

%% e: pearson correlation between same cells in diff track

load('D:\final_HDAC_AD_automatic_processing\dataloc\Remapping_linear_track_053019_053119_merged.mat')
load('Y:\Lujia\sfn2019\SFN2019 fig and text\Round17\figure3_lt_group.mat');
foldername=unique(destination);

behav_list={
'Horizontal1_Behav.mat';
'Horizontal2_Behav.mat';
'Vertical_Behav.mat';
};

for tk=1:6
   rateMap_corr1=rateMap_correlation(behavpos_half_ratemap_lt{tk,1}{1},behavpos_half_ratemap_lt{tk,1}{2},behavpos_half_ct_lt{tk,1}{1},behavpos_half_ct_lt{tk,1}{2},1);
   rateMap_corr2=rateMap_correlation(behavpos_half_ratemap_lt{tk,2}{1},behavpos_half_ratemap_lt{tk,2}{2},behavpos_half_ct_lt{tk,2}{1},behavpos_half_ct_lt{tk,2}{2},1);
   rateMap_corr3=rateMap_correlation(behavpos_half_ratemap_lt{tk,3}{1},behavpos_half_ratemap_lt{tk,3}{2},behavpos_half_ct_lt{tk,3}{1},behavpos_half_ct_lt{tk,3}{2},1);

   rateMap_corr_all_halves_lt(tk,1)=mean([nanmean(rateMap_corr1),nanmean(rateMap_corr2),nanmean(rateMap_corr3)]);
   
    parfor i=1:1000
        idx_rand_fr1=randperm(length(behavpos_half_ratemap_lt{tk,1}{1}));
        fr1_shuf=behavpos_half_ratemap_lt{tk,1}{1}(idx_rand_fr1);
        ct1_shuf=behavpos_half_ct_lt{tk,1}{1};
        
        idx_rand_fr2=randperm(length(behavpos_half_ratemap_lt{tk,1}{2}));
        fr2_shuf=behavpos_half_ratemap_lt{tk,1}{2}(idx_rand_fr2);
        ct2_shuf=behavpos_half_ct_lt{tk,1}{2};
        
        idx_rand_fr3=randperm(length(behavpos_half_ratemap_lt{tk,2}{1}));
        fr3_shuf=behavpos_half_ratemap_lt{tk,2}{1}(idx_rand_fr3);
        ct3_shuf=behavpos_half_ct_lt{tk,2}{1};
        
        idx_rand_fr4=randperm(length(behavpos_half_ratemap_lt{tk,2}{2}));
        fr4_shuf=behavpos_half_ratemap_lt{tk,2}{2}(idx_rand_fr4);
        ct4_shuf=behavpos_half_ct_lt{tk,2}{2};
        
        idx_rand_fr5=randperm(length(behavpos_half_ratemap_lt{tk,3}{1}));
        fr5_shuf=behavpos_half_ratemap_lt{tk,3}{1}(idx_rand_fr5);
        ct5_shuf=behavpos_half_ct_lt{tk,3}{1};
        
        idx_rand_fr6=randperm(length(behavpos_half_ratemap_lt{tk,3}{2}));
        fr6_shuf=behavpos_half_ratemap_lt{tk,3}{2}(idx_rand_fr6);
        ct6_shuf=behavpos_half_ct_lt{tk,3}{2};
        
       rateMap_corr1_shuf=rateMap_correlation(fr1_shuf,fr2_shuf,ct1_shuf,ct2_shuf,1);
       rateMap_corr2_shuf=rateMap_correlation(fr3_shuf,fr4_shuf,ct3_shuf,ct4_shuf,1);
       rateMap_corr3_shuf=rateMap_correlation(fr5_shuf,fr6_shuf,ct5_shuf,ct6_shuf,1);
       rateMap_corr_all_shuf_halves_lt(tk,1,i)=mean([nanmean(rateMap_corr1_shuf),nanmean(rateMap_corr2_shuf),nanmean(rateMap_corr3_shuf)]);
    end
toc;
end
rateMap_corr_all_shuf_halves_lt_1=quantile(rateMap_corr_all_shuf_halves_lt,.95,3);


load('D:\final_HDAC_AD_automatic_processing\dataloc\Remapping_linear_track_053019_053119_merged.mat')
foldername=unique(destination);
rateMap_corr_all=[];
rateMap_corr_all_shuf=[];
tic; 
for tk=1:6
    load([foldername{tk},'\neuronIndividuals_new.mat']);
    load([foldername{tk},'\Horizontal1_Behav.mat']);
    load([foldername{tk},'\thresh_and_ROI.mat'])
    [nnew1,bp_new,bt_new]=linearTrack_remove_water(neuronIndividuals_new{1},behav.position,behav.time);
    behavpos_all{1}=bp_new;
    behavtime_all{1}=bt_new;    
    load([foldername{tk},'\Horizontal2_Behav.mat']);
    [nnew2,bp_new,bt_new]=linearTrack_remove_water(neuronIndividuals_new{2},behav.position,behav.time);
    behavpos_all{2}=bp_new;
    behavtime_all{2}=bt_new;
    load([foldername{tk},'\Vertical_Behav.mat']);
    [nnew3,bp_new,bt_new]=linearTrack_remove_water(neuronIndividuals_new{3},behav.position,behav.time);
    behavpos_all{3}=bp_new;
    behavtime_all{3}=bt_new;
    
    [fr1,fr2,fr3,ct1,ct2,ct3]=figure3_loadFiringRate_C_recompute({nnew1,nnew2,nnew3},behavpos_all,behavtime_all,maxbehavROI);
    
    rateMap_corr1=rateMap_correlation(fr1,fr2,ct1,ct2,1);
    rateMap_corr2=rateMap_correlation(fr2,fr3,ct2,ct3,1);  
    rateMap_corr3=rateMap_correlation(fr1,fr3,ct1,ct3,1);  
    rateMap_corr_all(tk,:)=[nanmean(rateMap_corr1),nanmean(rateMap_corr2),nanmean(rateMap_corr3)];
   
    parfor i=1:1000
        idx_rand_fr1=randperm(length(fr1));
        fr1_shuf=fr1(idx_rand_fr1);
        ct1_shuf=ct1;
        
        idx_rand_fr2=randperm(length(fr2));
        fr2_shuf=fr2(idx_rand_fr2);
        ct2_shuf=ct2;
        
        idx_rand_fr3=randperm(length(fr3));
        fr3_shuf=fr3(idx_rand_fr3);
        ct3_shuf=ct3;
        
        rateMap_corr1_shuf=rateMap_correlation(fr1_shuf,fr2_shuf,ct1_shuf,ct2_shuf,1);
        rateMap_corr2_shuf=rateMap_correlation(fr2_shuf,fr3_shuf,ct2_shuf,ct3_shuf,1);  
        rateMap_corr3_shuf=rateMap_correlation(fr1_shuf,fr3_shuf,ct1_shuf,ct3_shuf,1);  
        rateMap_corr_all_shuf(tk,:,i)=[nanmean(rateMap_corr1_shuf),nanmean(rateMap_corr2_shuf),nanmean(rateMap_corr3_shuf)];
    end
toc;
end
rateMap_corr_all=[rateMap_corr_all_halves_lt,rateMap_corr_all];
rateMap_corr_all_shuf1=[rateMap_corr_all_shuf_halves_lt_1,quantile(rateMap_corr_all_shuf,.95,3)];

%% d: cluster overlap linear track
shuffle=1000;
h1_h2_adjrand=[];
h2_v_adjrand=[];
h1_v_adjrand=[];
h1_h2_shuffle_adjrand=[];
h2_v_shuffle_adjrand=[];
h1_v_shuffle_adjrand=[];

h1_h2_rand=[];
h2_v_rand=[];
h1_v_rand=[];
h1_h2_shuffle_rand=[];
h2_v_shuffle_rand=[];
h1_v_shuffle_rand=[];

clust_idx=4;
mice_all=[1:6]
for i=1:length(mice_all)
    mice=mice_all(i);
    gp=group(mice,:);
    [h1_h2_rand(i,1)]=new_cluster_overlap(gp{1}{clust_idx},gp{2}{clust_idx});
    [h2_v_rand(i,1)]=new_cluster_overlap(gp{2}{clust_idx},gp{3}{clust_idx});
    [h1_v_rand(i,1)]=new_cluster_overlap(gp{1}{clust_idx},gp{3}{clust_idx});

    for j=1:shuffle
        g1=gp{1}{clust_idx}(randperm(length(gp{1}{clust_idx})));
        g2=gp{2}{clust_idx}(randperm(length(gp{2}{clust_idx})));
        g3=gp{3}{clust_idx}(randperm(length(gp{3}{clust_idx})));
        [h1_h2_shuffle_rand(i,2)]=new_cluster_overlap(g1,g2);
        [h2_v_shuffle_rand(i,2)]=new_cluster_overlap(g2,g3);
        [h1_v_shuffle_rand(i,2)]=new_cluster_overlap(g1,g3);
    end
end

h1_h1_rand_half=[];
h1_h1_shuffle_rand_half=[];
clust_idx=4;
mice_all=[1:6];
for i=1:length(mice_all)
    mice=mice_all(i);
    group_half=group_record_half_lt(mice,:);
    [h1_h1_rand_half(i,1)]=new_cluster_overlap(group_half{1}{clust_idx}{1},group_half{1}{clust_idx}{2});

    for j=1:1000
        g1=group_half{1}{clust_idx}{1}(randperm(length(group_half{1}{clust_idx}{1})));
        g2=group_half{1}{clust_idx}{2}(randperm(length(group_half{1}{clust_idx}{2})));
        [h1_h1_shuffle_rand_half(i,2)]=new_cluster_overlap(g1,g2);
    end
end
dat_shuf=[quantile(h1_h1_shuffle_rand_half,0.95,2),quantile(h1_h2_shuffle_rand,0.95,2),quantile(h2_v_shuffle_rand,0.95,2),quantile(h1_v_shuffle_rand,0.95,2)];
dat_1=[h1_h1_rand_half,h1_h2_rand,h2_v_rand,h1_v_rand];
% quantile(dat(:),.95)




%% appendix: cluster calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% linear track first hor temporal half
load('D:\final_HDAC_AD_automatic_processing\dataloc\Remapping_linear_track_053019_053119_merged.mat')

foldername=unique(destination);
behavpos_half_ratemap_lt={};
behavpos_half_ct_lt={};
group_record_half_lt={};

behav_list={
'Horizontal1_Behav.mat';
'Horizontal2_Behav.mat';
'Vertical_Behav.mat';
};

for tk=1:6
    for j=1:3 
        for k=2:10
            load([foldername{tk},'\','neuronIndividuals_new.mat']);
            load([foldername{tk},'\','thresh_and_ROI.mat']);
            load([foldername{tk},'\',behav_list{j}]);
            [nnew,behavpos_half{tk,j},behavtime_half{tk,j}]=linearTrack_remove_water_cut_half_temporal(neuronIndividuals_new{j},behav.position,behav.time);
            [firingRateAll1,~,~,countTime1] = calculatingCellSpatialForSingleData_Suoqin(nnew{1},behavpos_half{tk,j}{1},behavtime_half{tk,j}{1},maxbehavROI,10,1:size(nnew{1}.C,1),3*std(nnew{1}.S,[],2),'S',[],[],[0 inf]);
            [firingRateAll2,~,~,countTime2] = calculatingCellSpatialForSingleData_Suoqin(nnew{2},behavpos_half{tk,j}{2},behavtime_half{tk,j}{2},maxbehavROI,10,1:size(nnew{2}.C,1),3*std(nnew{2}.S,[],2),'S',[],[],[0 inf]);
            behavpos_half_ratemap_lt{tk,j}={firingRateAll1,firingRateAll2};
            behavpos_half_ct_lt{tk,j}={countTime1,countTime2};
            group={};
            for h=1:2
                [group{h},CM,Z,cophList,cophList_adjusted]=cluster_determine_by_suoqin_NMF_firstPeakCoph_stability(nnew{h},100,10,k);
            end
            group_record_half_lt{tk,j}{k}=group;
        end 
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% linear track full
group={};
behav_list={
'Horizontal1_Behav.mat';
'Horizontal2_Behav.mat';
'Vertical_Behav.mat';
};
for tk=1:6
    for j=1:3
        for k=2:10
            load([foldername{tk},'\','neuronIndividuals_new.mat']);
            load([foldername{tk},'\',behav_list{j}]);
            [nnew]=linearTrack_remove_water(neuronIndividuals_new{j},behav.position,behav.time);
            [group,CM,Z,cophList,cophList_adjusted]=cluster_determine_by_suoqin_NMF_firstPeakCoph_stability(nnew,100,10,k);

            group{j,tk}{k}=group;
        end 
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% OF temporal half
load('D:\final_HDAC_AD_automatic_processing\dataloc\Remapping_square_circle_triangle_061119_061319_merged.mat')

foldername=unique(destination);
behavpos_half_ratemap_of={};
behavpos_half_ct_of={};
group_record_half_of={};

behav_list={
'triangle1_Behav.mat';
'circle1_Behav.mat';
'square1_Behav.mat';
'circle2_Behav.mat';
'square2_Behav.mat';
'triangle2_Behav.mat';
};

for tk=1:6
    for j=1:6 
        for k=2:10
            load([foldername{tk},'\','neuronIndividuals_new.mat']);
            load([foldername{tk},'\','thresh_and_ROI.mat']);
            load([foldername{tk},'\',behav_list{j}]);
            [nnew,behavpos_half{tk,j},behavtime_half{tk,j}]=OF_cut_half_temporal(neuronIndividuals_new{j},behav.position,behav.time);
            [firingRateAll1,~,~,countTime1] = calculatingCellSpatialForSingleData_Suoqin(nnew{1},behavpos_half{tk,j}{1},behavtime_half{tk,j}{1},maxbehavROI,10,1:size(nnew{1}.C,1),3*std(nnew{1}.S,[],2),'S',[],[],[0 inf]);
            [firingRateAll2,~,~,countTime2] = calculatingCellSpatialForSingleData_Suoqin(nnew{2},behavpos_half{tk,j}{2},behavtime_half{tk,j}{2},maxbehavROI,10,1:size(nnew{2}.C,1),3*std(nnew{2}.S,[],2),'S',[],[],[0 inf]);
            behavpos_half_ratemap_of{tk,j}={firingRateAll1,firingRateAll2};
            behavpos_half_ct_of{tk,j}={countTime1,countTime2};
            group={};
            for h=1:2
                [group{h},CM,Z,cophList,cophList_adjusted]=cluster_determine_by_suoqin_NMF_firstPeakCoph_stability(nnew{h},100,10,k);
            end
            group_record_half_of{tk,j}{k}=group;
        end 
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% OF full
group={};
behav_list={
'triangle1_Behav.mat';
'circle1_Behav.mat';
'square1_Behav.mat';
'circle2_Behav.mat';
'square2_Behav.mat';
'triangle2_Behav.mat';
};

for tk=5:6
    for j=1:6
        for k=2:10
            load([foldername{tk},'\','neuronIndividuals_new.mat']);
            load([foldername{tk},'\',behav_list{j}]);
            [group,CM,Z,cophList,cophList_adjusted]=cluster_determine_by_suoqin_NMF_firstPeakCoph_stability(neuronIndividuals_new{j},100,10,k);

            group{j,tk}{k}=group;
        end 
    end
end

%% appendix 2: manually assign color for fig3a
figure;
mice_all=[2]; %M3412
load('C:\Users\exx\Desktop\HDAC paper fig and method\SFN2019\SFN2019 fig\SCI clustNum plot\tri cir sqr\tri_cir_sqr_dataset.mat');

mice=mice_all(1)
load([foldername{mice},'\','neuronIndividuals_new.mat'])

for j=1:length(neuronIndividuals_new)
    [group{j}]=cluster_determine_by_suoqin_NMF_firstPeakCoph_stability(neuronIndividuals_new{j},100,10,6);
end

load('Y:\Lujia\sfn2019\SFN2019 fig and text\Round16_manuscriptXu0726\fig3 additional panels\panel a footprint manual color.mat');
gp_mapping=[[1:6];[1 2 6 5 4 3];[1:6];[1 5 3 4 6 2];[1 6 5 3 2 4];[1 4 3 2 5 6]];
A_color={};
A_color_region={};
for j=1:length(neuronIndividuals_new)
    gpt=group{j};
    gpt1=gpt;
    gpt(gpt1==1)=gp_mapping(j,1);
    gpt(gpt1==2)=gp_mapping(j,2);
    gpt(gpt1==3)=gp_mapping(j,3);
    gpt(gpt1==4)=gp_mapping(j,4);
    gpt(gpt1==5)=gp_mapping(j,5);
    gpt(gpt1==6)=gp_mapping(j,6);
    [A_color{j},A_color_region{j}]=DBSCAN_region_quantify_func_single_group_no_plot(gpt,neuronIndividuals_new);

    subplot(2,6,j)
    imagesc(A_color{j});
    subplot(2,6,j+6)
    imagesc(A_color_region{j});    
end
