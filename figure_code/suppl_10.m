%% suppl 13: ICA - kmean based consensus compare
color_clust=distinguishable_colors(20);
run('D:\Xu_clusterting_paper_prep11_2020\final_code\data_prepare\neuron_data_info.m')
all_colors=distinguishable_colors(10);
load(['D:\Xu_clusterting_paper_prep11_2020\final_code\final_cluster_data\cluster_optimal_num\Fig2_cir_rec_clust_original.mat'])


%% part 1AI163 fig3

% k-mean clustering
tic;
for j=1:5
    for j1=1:7
        load([foldername_AI163{j},'\','neuronIndividuals_new.mat'])
        [~,group_ori_AI163{j,j1}]=cluster_determine_by_suoqin_NMF_firstPeakCoph_022422(neuronIndividuals_new{j1},100,10,[]);
        [d1,d2]=size(neuronIndividuals_new{1}.Cn);
        A_color_ori_AI163{j,j1}=cluster_spatial_footprint_colormap(neuronIndividuals_new,d1,d2,color_clust,group_ori_AI163{j,j1},0.75);
        disp('finish')
        toc;
    end
end

% Torts 2013 ICA
color_clust=distinguishable_colors(20);

for j=1:5
    for j1=1:7
        load([foldername_AI163{j},'\','neuronIndividuals_new.mat'])
        AssemblyTemplates_AI163{j,j1} = assembly_patterns_num(neuronIndividuals_new{j1}.C,max(group_ori_AI163{j,j1}));
        [time_projection_AI163{j,j1}] = assembly_activity(AssemblyTemplates_AI163{j,j1},neuronIndividuals_new{j1}.C);
        
        group_tort_AI163{j,j1}=AssemblyTemplateCellGroupInfer(AssemblyTemplates_AI163{j,j1},time_projection_AI163{j,j1},neuronIndividuals_new{j1}.C,[]);
        [~,group_tort_AI163{j,j1},~]=alignClusterIdx(group_ori_AI163{j,j1},group_tort_AI163{j,j1});
        
        [d1,d2]=size(neuronIndividuals_new{1}.Cn);
        A_color_tort_AI163{j,j1}=cluster_spatial_footprint_colormap(neuronIndividuals_new,d1,d2,color_clust,group_tort_AI163{j,j1},0.75);
        
%         group_tort_high_AI163{j,j1}=AssemblyTemplateCellGroupInfer(AssemblyTemplates_AI163{j,j1},time_projection_AI163{j,j1},neuronIndividuals_new{j1}.C,0.90);
%         [~,group_tort_high_AI163{j,j1},~]=alignClusterIdx(group_ori_AI163{j,j1},group_tort_high_AI163{j,j1});
%         
%         [d1,d2]=size(neuronIndividuals_new{1}.Cn);
%         A_color_tort_high_AI163{j,j1}=cluster_spatial_footprint_colormap(neuronIndividuals_new,d1,d2,color_clust,group_tort_high_AI163{j,j1},0.75);

        disp('finish')
    end
end

toc;

%% part2 multigeo

% k-mean clustering
for j=1:6
    for j1=1:6
        load([foldername_multiGeo{j},'\','neuronIndividuals_new.mat'])
        [~,group_ori_multiGeo{j,j1}]=cluster_determine_by_suoqin_NMF_firstPeakCoph_022422(neuronIndividuals_new{j1},100,10,[]);
        [d1,d2]=size(neuronIndividuals_new{1}.Cn);
        A_color_ori_multiGeo{j,j1}=cluster_spatial_footprint_colormap(neuronIndividuals_new,d1,d2,color_clust,group_ori_multiGeo{j,j1},0.75);
        disp('finish')
    end
end

% Torts 2013 ICA

for j=1:6
    for j1=1:6
        load([foldername_multiGeo{j},'\','neuronIndividuals_new.mat'])
        AssemblyTemplates_multiGeo{j,j1} = assembly_patterns_num(neuronIndividuals_new{j1}.C,max(group_ori_multiGeo{j,j1}));
        [time_projection_multiGeo{j,j1}] = assembly_activity(AssemblyTemplates_multiGeo{j,j1},neuronIndividuals_new{j1}.C);
        
        group_tort_multiGeo{j,j1}=AssemblyTemplateCellGroupInfer(AssemblyTemplates_multiGeo{j,j1},time_projection_multiGeo{j,j1},neuronIndividuals_new{j1}.C,[]);  
        [~,group_tort_multiGeo{j,j1},~]=alignClusterIdx(group_ori_multiGeo{j,j1},group_tort_multiGeo{j,j1});
        
        [d1,d2]=size(neuronIndividuals_new{1}.Cn);
        A_color_tort_multiGeo{j,j1}=cluster_spatial_footprint_colormap(neuronIndividuals_new,d1,d2,color_clust,group_tort_multiGeo{j,j1},0.75);
        
%         group_tort_high_multiGeo{j,j1}=AssemblyTemplateCellGroupInfer(AssemblyTemplates_multiGeo{j,j1},time_projection_multiGeo{j,j1},neuronIndividuals_new{j1}.C,0.50);
%         [~,group_tort_high_multiGeo{j,j1},~]=alignClusterIdx(group_ori_multiGeo{j,j1},group_tort_high_multiGeo{j,j1});
%         
%         [d1,d2]=size(neuronIndividuals_new{1}.Cn);
%         A_color_tort_high_multiGeo{j,j1}=cluster_spatial_footprint_colormap(neuronIndividuals_new,d1,d2,color_clust,group_tort_high_multiGeo{j,j1},0.75);
% 
        disp('finish')
    end
end

%% part3 fig2

% k-mean clustering
for j=1:12
    for j1=1:2
        load([foldername_fig2{j},'\','neuronIndividuals_new.mat'])
        [~,group_ori_fig2{j,j1}]=cluster_determine_by_suoqin_NMF_firstPeakCoph_022422(neuronIndividuals_new{j1},100,10,7);
        [d1,d2]=size(neuronIndividuals_new{1}.Cn);
        A_color_ori_fig2{j,j1}=cluster_spatial_footprint_colormap(neuronIndividuals_new,d1,d2,color_clust,group_ori_fig2{j,j1},0.75);
        disp('finish')
    end
end

% Torts 2013 ICA
for j=1:12
    for j1=1:2
        load([foldername_fig2{j},'\','neuronIndividuals_new.mat'])
        AssemblyTemplates_fig2{j,j1} = assembly_patterns_num(neuronIndividuals_new{j1}.C,max(group_ori_fig2{j,j1}));
        [time_projection_fig2{j,j1}] = assembly_activity(AssemblyTemplates_fig2{j,j1},neuronIndividuals_new{j1}.C);
        
        group_tort_fig2{j,j1}=AssemblyTemplateCellGroupInfer(AssemblyTemplates_fig2{j,j1},time_projection_fig2{j,j1},neuronIndividuals_new{j1}.C,[]);
        [~,group_tort_fig2{j,j1},~]=alignClusterIdx(group_ori_fig2{j,j1},group_tort_fig2{j,j1});
        
        [d1,d2]=size(neuronIndividuals_new{1}.Cn);
        A_color_tort_fig2{j,j1}=cluster_spatial_footprint_colormap(neuronIndividuals_new,d1,d2,color_clust,group_tort_fig2{j,j1},0.75);
        
%         group_tort_high_fig2{j,j1}=AssemblyTemplateCellGroupInfer(AssemblyTemplates_fig2{j,j1},time_projection_fig2{j,j1},neuronIndividuals_new{j1}.C,0.90);
%         [~,group_tort_high_fig2{j,j1},~]=alignClusterIdx(group_ori_fig2{j,j1},group_tort_high_fig2{j,j1});
%         
%         [d1,d2]=size(neuronIndividuals_new{1}.Cn);
%         A_color_tort_high_fig2{j,j1}=cluster_spatial_footprint_colormap(neuronIndividuals_new,d1,d2,color_clust,group_tort_high_fig2{j,j1},0.75);
% 
    end
end

%% illustration
% part 1: time series
% multiGeo:
load([foldername_multiGeo{6},'\','neuronIndividuals_new.mat']);
clustered_timeSeries2(neuronIndividuals_new{3}.C,group_ori_multiGeo{6,3});
clustered_timeSeries2(neuronIndividuals_new{3}.C,group_tort_multiGeo{6,3});


% AI163
load([foldername_AI163{1},'\','neuronIndividuals_new.mat']);
clustered_timeSeries2(neuronIndividuals_new{1}.C,group_ori_AI163{1,1});
clustered_timeSeries2(neuronIndividuals_new{1}.C,group_tort_AI163{1,1});

% Fig2
load([foldername_fig2{1},'\','neuronIndividuals_new.mat']);
clustered_timeSeries2(neuronIndividuals_new{1}.C,group_ori_fig2{1,1});
clustered_timeSeries2(neuronIndividuals_new{1}.C,group_tort_fig2{1,1});

% part 2: assembly series
% multiGeo:
load([foldername_multiGeo{6},'\','neuronIndividuals_new.mat']);
assembly_series_illustrate(neuronIndividuals_new{3}.C,group_ori_multiGeo{6,3});
assembly_series_illustrate(neuronIndividuals_new{3}.C,group_tort_multiGeo{6,3});

% AI163
load([foldername_AI163{1},'\','neuronIndividuals_new.mat']);
assembly_series_illustrate(neuronIndividuals_new{1}.C,group_ori_AI163{1,1});
assembly_series_illustrate(neuronIndividuals_new{1}.C,group_tort_AI163{1,1});

% Fig2
load([foldername_fig2{1},'\','neuronIndividuals_new.mat']);
assembly_series_illustrate(neuronIndividuals_new{1}.C,group_ori_fig2{1,1});
assembly_series_illustrate(neuronIndividuals_new{1}.C,group_tort_fig2{1,1});

% part 3: foot prints
subplot(121);imagesc(A_color_ori_multiGeo{6,3});subplot(122);imagesc(A_color_tort_multiGeo{6,3})
subplot(121);imagesc(A_color_ori_AI163{1,1});subplot(122);imagesc(A_color_tort_AI163{1,1})
subplot(121);imagesc(A_color_ori_fig2{1,1});subplot(122);imagesc(A_color_tort_fig2{1,1})

%% overlap analysis

all_clust_overlap_ai163=[];
for i=1:size(group_ori_AI163,1)
    for j=1:size(group_ori_AI163,2)
        g1=group_ori_AI163{i,j};
        g2=group_tort_AI163{i,j};
        del_idx=(double(g1==-1)+double(g2==-1))>0;
        g1(del_idx)=[];
        g2(del_idx)=[];
        all_clust_overlap_ai163(i,j)=new_cluster_overlap_latest(g1,g2);
    end
end

all_clust_overlap_multiGeo=[];
for i=1:size(group_ori_multiGeo,1)
    for j=1:size(group_ori_multiGeo,2)
        g1=group_ori_multiGeo{i,j};
        g2=group_tort_multiGeo{i,j};
        del_idx=(double(g1==-1)+double(g2==-1))>0;
        g1(del_idx)=[];
        g2(del_idx)=[];
        all_clust_overlap_multiGeo(i,j)=new_cluster_overlap_latest(g1,g2);
    end
end

all_clust_overlap_fig2=[];
for i=1:size(group_ori_fig2,1)
    for j=1:size(group_ori_fig2,2)
        g1=group_ori_fig2{i,j};
        g2=group_tort_fig2{i,j};
        del_idx=(double(g1==-1)+double(g2==-1))>0;
        g1(del_idx)=[];
        g2(del_idx)=[];
        all_clust_overlap_fig2(i,j)=new_cluster_overlap_latest(g1,g2);
    end
end

%% intra_cluster correlation, inter_cluster correlation, and compare with k_mean
for tk=1:length(foldername_multiGeo)
    load([foldername_multiGeo{tk},'\neuronIndividuals_new.mat']);
    [~,~,~,~,~,~,~,intra_all_multiGeo{tk},inter_all_multiGeo{tk},intra_shuffle_all_multiGeo{tk}]=intra_inter_cluster_corr_dis(neuronIndividuals_new,group_ori_multiGeo{tk,1},1,'corr');
    [~,~,~,~,~,~,~,intra_all_multiGeo_tort{tk},inter_all_multiGeo_tort{tk},intra_shuffle_all_multiGeo_tort{tk}]=intra_inter_cluster_corr_dis(neuronIndividuals_new,group_tort_multiGeo{tk,1},1,'corr');

end

for tk=1:length(foldername_fig2)
    load([foldername_fig2{tk},'\neuronIndividuals_new.mat']);
    [~,~,~,~,~,~,~,intra_all_fig2{tk},inter_all_fig2{tk},intra_shuffle_all_fig2{tk}]=intra_inter_cluster_corr_dis(neuronIndividuals_new,group_ori_fig2{tk,1},1,'corr');
    [~,~,~,~,~,~,~,intra_all_fig2_tort{tk},inter_all_fig2_tort{tk},intra_shuffle_all_fig2_tort{tk}]=intra_inter_cluster_corr_dis(neuronIndividuals_new,group_tort_fig2{tk,1},1,'corr');

end

for tk=1:length(foldername_AI163)
    load([foldername_AI163{tk},'\neuronIndividuals_new.mat']);
    [~,~,~,~,~,~,~,intra_all_AI163{tk},inter_all_AI163{tk},intra_shuffle_all_AI163{tk}]=intra_inter_cluster_corr_dis(neuronIndividuals_new,group_ori_AI163{tk,1},1,'corr');
    [~,~,~,~,~,~,~,intra_all_AI163_tort{tk},inter_all_AI163_tort{tk},intra_shuffle_all_AI163_tort{tk}]=intra_inter_cluster_corr_dis(neuronIndividuals_new,group_tort_AI163{tk,1},1,'corr');
end

colorss=colormap(lines);
colorss=[1,0,1;colorss(4:6,:)];
figure;
subplot(131); hold on
cluster_correlation_cdf_plot({intra_all_multiGeo,intra_all_multiGeo_tort,inter_all_multiGeo_tort,intra_shuffle_all_multiGeo_tort},colorss);
subplot(132); hold on
cluster_correlation_cdf_plot({intra_all_fig2,intra_all_fig2_tort,inter_all_fig2_tort,intra_shuffle_all_fig2_tort},colorss);
subplot(133); hold on
cluster_correlation_cdf_plot({intra_all_AI163,intra_all_AI163_tort,inter_all_AI163_tort,intra_shuffle_all_AI163_tort},colorss);

% statistic
[v1,v2,v3,v4,p11,p21,p31]=suppl12_intra_inter_statistic(intra_all_multiGeo,intra_all_multiGeo_tort,intra_shuffle_all_multiGeo_tort,inter_all_multiGeo_tort);
[v1,v2,v3,v4,p12,p22,p32]=suppl12_intra_inter_statistic(intra_all_fig2,intra_all_fig2_tort,intra_shuffle_all_fig2_tort,inter_all_fig2_tort);
[v1,v2,v3,v4,p13,p23,p33]=suppl12_intra_inter_statistic(intra_all_AI163,intra_all_AI163_tort,intra_shuffle_all_AI163_tort,inter_all_AI163_tort);


%% spatial cluster region size
% change pattern, y axis
[avg_region_kmean_multiGeo,avg_region_shuf_kmean_multiGeo,avg_region_ICA_multiGeo,avg_region_ICA_shuf_multiGeo]=clustNum_ICA_kmean(foldername_multiGeo,3);
[avg_region_kmean_fig2,avg_region_shuf_kmean_fig2,avg_region_ICA_fig2,avg_region_ICA_shuf_fig2]=clustNum_ICA_kmean(foldername_fig2,1);
[avg_region_kmean_AI163,avg_region_shuf_kmean_AI163,avg_region_ICA_AI163,avg_region_ICA_shuf_AI163]=clustNum_ICA_kmean(foldername_AI163,1);

subplot(131)
plot_cluster_region_with_shuf(avg_region_ICA_multiGeo(:,3),avg_region_ICA_shuf_multiGeo(:,3),[0,9000]);
subplot(132)
plot_cluster_region_with_shuf(avg_region_ICA_fig2(:,1),avg_region_ICA_shuf_fig2(:,1),[0,4500]);
subplot(133)
plot_cluster_region_with_shuf(avg_region_ICA_AI163(:,1),avg_region_ICA_shuf_AI163(:,1),[0,2500]);

%% intra-cluster-inter cluster correlation over distance
[all_intra_mGeo,all_inter_mGeo,~,~,f_mGeo,gof_mGeo]=pairwise_dis_tempCorr_031722(foldername_multiGeo,group_tort_multiGeo,3)
[all_intra_fig2,all_inter_fig2,~,~,f_fig2,gof_fig2]=pairwise_dis_tempCorr_031722(foldername_fig2,group_tort_fig2,1)
[all_intra_ai163,all_inter_ai163,~,~,f_ai163,gof_ai163]=pairwise_dis_tempCorr_031722(foldername_AI163,group_tort_AI163,1)

[H, pValue1, KSstatistic] = kstest_2s_2d(all_intra_mGeo(1:end,:), all_inter_mGeo(1:end,:));
[H, pValue2, KSstatistic] = kstest_2s_2d(all_intra_fig2(1:end,:), all_inter_fig2(1:end,:));
[H, pValue3, KSstatistic] = kstest_2s_2d(all_intra_ai163(1:end,:), all_inter_ai163(1:end,:));

all_x=[all_intra_mGeo(:,1);all_inter_mGeo(:,1)];
x=linspace( min(all_x), max(all_x), 100)';
y_intra=f_mGeo{1}.a*x.^f_mGeo{1}.b+f_mGeo{1}.c;
y_inter=f_mGeo{2}.a*x.^f_mGeo{2}.b+f_mGeo{2}.c;
signrank(y_intra,y_inter)

all_x=[all_intra_fig2(:,1);all_inter_fig2(:,1)];
x=linspace( min(all_x), max(all_x), 100)';
y_intra=f_fig2{1}.a*x.^f_fig2{1}.b+f_fig2{1}.c;
y_inter=f_fig2{2}.a*x.^f_fig2{2}.b+f_fig2{2}.c;
signrank(y_intra,y_inter)

all_x=[all_intra_ai163(:,1);all_inter_ai163(:,1)];
x=linspace( min(all_x), max(all_x), 100)';
y_intra=f_ai163{1}.a*x.^f_ai163{1}.b+f_ai163{1}.c;
y_inter=f_ai163{2}.a*x.^f_ai163{2}.b+f_ai163{2}.c;
signrank(y_intra,y_inter)

disp([mean(all_intra_mGeo(1:end,1)), sem(all_intra_mGeo(1:end,1),1),mean(all_inter_mGeo(1:end,1)), sem(all_inter_mGeo(1:end,1),1), pValue1]);
disp([mean(all_intra_fig2(1:end,1)), sem(all_intra_fig2(1:end,1),1),mean(all_inter_fig2(1:end,1)), sem(all_inter_fig2(1:end,1),1), pValue2]);
disp([mean(all_intra_ai163(1:end,1)), sem(all_intra_ai163(1:end,1),1),mean(all_inter_ai163(1:end,1)), sem(all_inter_ai163(1:end,1),1), pValue3]);

%% quantify the amount of neurons assigned as grey, in ICA
unassigned_fraction={[],[],[]};
for i=1:size(group_tort_multiGeo,1)
    unassigned_fraction{1}(i,1)=sum(group_tort_multiGeo{i,3}==-1)/length(group_tort_multiGeo{i,3});
end
for i=1:size(group_tort_fig2,1)
    unassigned_fraction{2}(i,1)=sum(group_tort_fig2{i,1}==-1)/length(group_tort_fig2{i,1});
end
for i=1:size(group_tort_AI163,1)
    unassigned_fraction{3}(i,1)=sum(group_tort_AI163{i,6}==-1)/length(group_tort_AI163{i,6});
end