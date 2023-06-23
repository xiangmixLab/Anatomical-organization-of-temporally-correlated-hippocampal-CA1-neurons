%% supplemental Figure 4: PC distribution
load(['D:\Xu_clusterting_paper_prep11_2020\final_code\final_cluster_data\cluster_optimal_num\Fig1_3_multiGeo_clust_original.mat'])

%% 1.PC calculation
binsize=10;

fldname={
    'D:\Xu_clusterting_paper_prep11_2020\arranged_final_data\Fig_1_3_tri_cir_sqr\M3411'		
    'D:\Xu_clusterting_paper_prep11_2020\arranged_final_data\Fig_1_3_tri_cir_sqr\M3412'		
    'D:\Xu_clusterting_paper_prep11_2020\arranged_final_data\Fig_1_3_tri_cir_sqr\M3421F'		
    'D:\Xu_clusterting_paper_prep11_2020\arranged_final_data\Fig_1_3_tri_cir_sqr\M3422F'		
    'D:\Xu_clusterting_paper_prep11_2020\arranged_final_data\Fig_1_3_tri_cir_sqr\M3424F'		
    'D:\Xu_clusterting_paper_prep11_2020\arranged_final_data\Fig_1_3_tri_cir_sqr\M3425F'		
};

nnew_all={};
all_behav={};
for i=1:length(fldname)
    load([fldname{i},'\','neuronIndividuals_new.mat']);
    load([fldname{i},'\','behav.mat']);
    for j=1:length(neuronIndividuals_new)
        all_behav{i,j}=behavIndividuals{j};
        all_behav{i,j}.VidObj=[];
        nnew_all{i}=neuronIndividuals_new;
    end
end

% PC
all_pc=cell(6,6);
all_infoscore=cell(6,6);
all_infoscore_norm=cell(6,6);
all_coherence=cell(6,6);
tic;
for i=1:6
    for j=2

        [place_cells,infoScore,infoScore_norm,coherencee] = permutingSpike_adapt_040821(nnew_all{i}{j},all_behav{i}.position,all_behav{i}.time,'S',0.1,10,5,'all',0.3);  

        all_pc{i,j}=place_cells;
        all_infoscore{i,j}=infoScore;
        all_infoscore_norm{i,j}=infoScore_norm;
        all_coherence{i,j}=coherencee;

    end
end
toc;

%% 2. illustration

colorClusters_all=distinguishable_colors(11);
colorClusters_all(11,:)=[0.75,0.75,0.75];
for j=1:length(nnew_all)
    for j1=2
        nnew_all{j}{j1}.imageSize=[240,376];
        nnew_all{j}{j1}.centroid=neuron_centroid_calculation(nnew_all{j}{j1},nnew_all{j}{j1}.imageSize);
        [A_color{j,j1},A_color_region{j,j1}]=DBSCAN_region_quantify_022422(group_ori_multiGeo{j,j1},nnew_all{j},[]);
        
        gp1=group_ori_multiGeo{j,j1};
        gp1(~ismember(1:length(gp1),all_pc{j,2}{2}))=11;
        A_color_pc{j,j1}=cluster_spatial_footprint_colormap({nnew_all{j}{j1}},240,376,colorClusters_all,gp1,0.7);

    end
end

midx=[1:6];

figure;
ctt=1;
for j=midx
    for j1=2
        subplot(6,2,j1+(ctt-1)*2);
        imagesc(A_color_pc{j,j1});        
    end
    ctt=ctt+1;
end

%% 3. pairwise distance

for tk=1:6
    [~,~,~,~,~,~,~,intra_all{tk},inter_all{tk},intra_shuffle_all{tk}]=intra_inter_cluster_corr_dis({nnew_all{tk}{2}},group_ori_multiGeo{j,2},1,'dis');
    
    pc_idx=all_pc{tk,2}{2};
    del_idx=find(~ismember(1:size(nnew_all{tk}{2}.C,1),pc_idx)==1);
    
    npc=nnew_all{tk}{2};
    npc.C(del_idx,:)=[];
    npc.C_raw(del_idx,:)=[];
    npc.S(del_idx,:)=[];
    npc.A(:,del_idx)=[];
    npc.centroid(del_idx,:)=[];
    gppc=group_ori_multiGeo{tk,2};
    gppc=gppc(pc_idx);
    
    [~,~,~,~,~,~,~,intra_all_pc{tk}]=intra_inter_cluster_corr_dis({npc},gppc,1,'dis');
    [~,~,~,~,~,~,~,intra_all_pc_all{tk}]=intra_inter_cluster_corr_dis({npc},double(gppc>0),1,'dis');
    [dis_all{tk}]=intra_inter_cluster_corr_dis({nnew_all{tk}{2}},double(group_ori_multiGeo{tk,2}>0),1,'dis');

end

% intra_all,intra_pc,all_pc
colorClusters_all=lines(256);

colorClusters_all_fig=[];
for i=1:256
    colorClusters_all_fig(1,i,:)=colorClusters_all(i,:);
end

figure;
h1=cdfplot(cell2mat(intra_all'));
hold on;
h2=cdfplot(cell2mat(intra_all_pc'));
h3=cdfplot(cell2mat(intra_all_pc_all'));
set(h1,'color',colorClusters_all(4,:));
set(h2,'color',colorClusters_all(5,:));
set(h3,'color',colorClusters_all(6,:));

p1=infer_cdf_loc(cell2mat(intra_all'),nanmean(cell2mat(intra_all')));
p2=infer_cdf_loc(cell2mat(intra_all_pc'),nanmean(cell2mat(intra_all_pc')));
p3=infer_cdf_loc(cell2mat(intra_all_pc_all'),nanmean(cell2mat(intra_all_pc_all')));

plot(p1(1),p1(2),'.','color',colorClusters_all(4,:),'MarkerSize',20)
plot(p2(1),p2(2),'.','color',colorClusters_all(5,:),'MarkerSize',20)
plot(p3(1),p3(2),'.','color',colorClusters_all(6,:),'MarkerSize',20)

k1=cell2mat(intra_all');
k2=cell2mat(intra_all_pc');
k3=cell2mat(intra_all_pc_all');


[H1,P1] = kstest2(k1,k2);
[H2,P2] = kstest2(k1,k3);
[H3,P3] = kstest2(k2,k3);

mean(cell2mat(intra_all')*2)% dis correction
sem(cell2mat(intra_all')*2,1)

mean(cell2mat(intra_all_pc')*2)
sem(cell2mat(intra_all_pc')*2,1)

mean(cell2mat(intra_all_pc_all')*2)
sem(cell2mat(intra_all_pc_all')*2,1)

mean(cell2mat(dis_all')*2)

%% 4. pairwise CORR
for tk=1:6
    [~,~,~,~,~,~,~,intra_all{tk},inter_all{tk},intra_shuffle_all{tk}]=intra_inter_cluster_corr_dis({nnew_all{tk}{2}},group_ori_multiGeo{tk,2},1,'corr');
    
    pc_idx=all_pc{tk,2}{2};
    del_idx=find(~ismember(1:size(nnew_all{tk}{2}.C,1),pc_idx)==1);
    
    npc=nnew_all{tk}{2};
    npc.C(del_idx,:)=[];
    npc.C_raw(del_idx,:)=[];
    npc.S(del_idx,:)=[];
    npc.A(:,del_idx)=[];
    npc.centroid(del_idx,:)=[];
    gppc=group_ori_multiGeo{tk,2};
    gppc=gppc(pc_idx);
    
    [~,~,~,~,~,~,~,intra_all_pc{tk}]=intra_inter_cluster_corr_dis({npc},gppc,1,'corr');
    [~,~,~,~,~,~,~,intra_all_pc_all{tk}]=intra_inter_cluster_corr_dis({npc},double(gppc>0),1,'corr');
    [corr_all{tk}]=intra_inter_cluster_corr_dis({nnew_all{tk}{2}},double(group_ori_multiGeo{tk,2}>0),1,'corr');
end

% intra_all,intra_pc,all_pc
colorClusters_all=lines(256);

colorClusters_all_fig=[];
for i=1:256
    colorClusters_all_fig(1,i,:)=colorClusters_all(i,:);
end
imagesc(colorClusters_all_fig);

figure;
h1=cdfplot(cell2mat(intra_all'));
hold on;
h2=cdfplot(cell2mat(intra_all_pc'));
h3=cdfplot(cell2mat(intra_all_pc_all'));
set(h1,'color',colorClusters_all(4,:));
set(h2,'color',colorClusters_all(5,:));
set(h3,'color',colorClusters_all(6,:));

p1=infer_cdf_loc(cell2mat(intra_all'),nanmean(cell2mat(intra_all')));
p2=infer_cdf_loc(cell2mat(intra_all_pc'),nanmean(cell2mat(intra_all_pc')));
p3=infer_cdf_loc(cell2mat(intra_all_pc_all'),nanmean(cell2mat(intra_all_pc_all')));

plot(p1(1),p1(2),'.','color',colorClusters_all(4,:),'MarkerSize',20)
plot(p2(1),p2(2),'.','color',colorClusters_all(5,:),'MarkerSize',20)
plot(p3(1),p3(2),'.','color',colorClusters_all(6,:),'MarkerSize',20)

[H1,P1] = kstest2(cell2mat(intra_all'),cell2mat(intra_all_pc'));
[H2,P2] = kstest2(cell2mat(intra_all'),cell2mat(intra_all_pc_all'));
[H3,P3] = kstest2(cell2mat(intra_all_pc'),cell2mat(intra_all_pc_all'));

k1=cell2mat(intra_all_pc');
k2=cell2mat(intra_all');
k3=cell2mat(intra_all_pc_all');

[H2,P2] = kstest2(k1(1:10:end),k3(1:10:end));

mean(k1)% dis correction
sem(k1,1)

mean(k2)
sem(k2,1)

mean(k3)
sem(k3,1)

