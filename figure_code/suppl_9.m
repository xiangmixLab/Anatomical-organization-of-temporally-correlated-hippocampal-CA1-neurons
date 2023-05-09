%% supplemental Fig. 8: Ai163 region analysis
run('D:\Xu_clusterting_paper_prep11_2020\final_code\data_prepare\neuron_data_info.m')

foldername=foldername_AI163;


%% illustrate behavior trace
load([foldername{3},'\behav.mat'])
plot(behavIndividuals{1}.position(:,1),behavIndividuals{1}.position(:,2),'lineWidth',2);

%% intra-inter cluster corr
load('D:\Xu_clusterting_paper_prep11_2020\final_code\final_cluster_data\cluster_optimal_num\Fig3_barrier_clust_original.mat')
for tk=1:length(foldername)
    load([foldername{tk},'\neuronIndividuals_new.mat']);
    [intra_all,inter_all,~,~,intra_shuffle_all,~,~,intra_all_neuron{tk},inter_all_neuron{tk},intra_shuffle_all_neuron{tk}]=intra_inter_cluster_corr_dis(neuronIndividuals_new,group_ori_AI163{tk,cond},cond,'corr');
end

colorss=colormap(lines);
colorss=[colorss(4:6,:)];
cluster_correlation_cdf_plot({{intra_all},{inter_all},{intra_shuffle_all}},colorss);
[v,p]=intra_inter_statistic_perNeuron({intra_all_neuron,inter_all_neuron,intra_shuffle_all_neuron},{'intra','inter','intra_shuffle'});

%% footprint
cond=8;

for i=1:5
    load([foldername{i},'\neuronIndividuals_new.mat']);
    [A_color1,A_color_region1]=DBSCAN_region_quantify_022422(group_ori_AI163{i,cond},neuronIndividuals_new,[]);
    subplot(2,5,i)
    imagesc(A_color1)
    subplot(2,5,i+5)
    imagesc(A_color_region1)
end
%% calcualte region, for differetn cluster num
load('D:\Xu_clusterting_paper_prep11_2020\final_code\final_cluster_data\cluster_region_sz\fig3_barrier_cluster_Region_sz.mat')

% current data is calculated as removing all neurons with max corr<=0.3
conds=[1,2,3,4,6,7,9]
avg_region_cond1_mat=[];
for i=1:size(avg_region_AI163_k,1)
    temp=[];
    for j=1:length(conds)
        temp=[temp;cell2mat(avg_region_AI163_k{i,conds(j)})];
    end
    avg_region_cond1_mat(i,:)=nanmean(temp,1); 
end

avg_region_cond1_shuf_mat=[];
for i=1:size(avg_region_shuf_AI163_k,1) 
    temp=[];
    for k=1:length(conds)
        shuf_region=avg_region_shuf_AI163_k{i,conds(k)};
        shuf_region_mat=[];
        for j=1:length(shuf_region)
            shuf_region_mat(j)=nanmean(shuf_region{j});
        end
        temp(k,:)=shuf_region_mat;
    end
    avg_region_cond1_shuf_mat(i,:)=nanmean(temp,1);
end

% illustrate single mice
for i=1:5
    subplot(1,5,i)
    plot(avg_region_cond1_mat(i,:));hold on;
    plot(avg_region_cond1_shuf_mat(i,2:end));
    signrank(avg_region_cond1_mat(i,:),avg_region_cond1_shuf_mat(i,2:end))
end

% illustrate all mice with error bar, ALL
errorbar(mean(avg_region_cond1_mat(:,:),1),std(avg_region_cond1_mat(:,:),[],1)/size(avg_region_cond1_mat,1).^0.5)
hold on;
errorbar(mean(avg_region_cond1_shuf_mat(:,2:end),1),std(avg_region_cond1_shuf_mat(:,2:end)/size(avg_region_cond1_mat,1).^0.5,[],1))


%% compare size with optimal cluster number for each animal
region_compare=[];
for i=1:size(group_ori_AI163,1)
    region_compare(i,:)=[avg_region_cond1_mat(i,max(group_ori_AI163{i,cond})-1),avg_region_cond1_shuf_mat(i,max(group_ori_AI163{i,cond}))];
end

ranksum(region_compare(:,1),region_compare(:,2))
signrank(region_compare(:,1),region_compare(:,2))

