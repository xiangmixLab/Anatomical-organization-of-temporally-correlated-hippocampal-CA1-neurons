%% AI163 footprint and size
foldername={
    'D:\Xu_clusterting_paper_prep11_2020\final_data\final_neuron_behav_data\Fig_3_barrier\101F\'		
    'D:\Xu_clusterting_paper_prep11_2020\final_data\final_neuron_behav_data\Fig_3_barrier\102F\'		
    'D:\Xu_clusterting_paper_prep11_2020\final_data\final_neuron_behav_data\Fig_3_barrier\103F\'		
    'D:\Xu_clusterting_paper_prep11_2020\final_data\final_neuron_behav_data\Fig_3_barrier\2833M\'		
    'D:\Xu_clusterting_paper_prep11_2020\final_data\final_neuron_behav_data\Fig_3_barrier\840F\'		
}


for j=1:length(foldername)
    load([foldername{j},'\','neuronIndividuals_new.mat'])
    for j1=7 % brick's empty
        [footprint{j,j1}{1},footprint{j,j1}{2},avg_region{j,j1},~,~,~,avg_region_shuffled{j,j1}]=DBSCAN_region_quantify_022422(group_ori_AI163{j,j1},neuronIndividuals_new,[]);
    end
end

%% plot

% ft
figure;
for i=1:5
    subplot(1,5,i);
    imagesc(footprint{i,7}{1});
end

% region
figure;
for i=1:5
    subplot(1,5,i);
    imagesc(footprint{i,7}{2});
end

%% size and shuffle size
load('D:\Xu_clusterting_paper_prep11_2020\final_data\final_cluster_data\diferent_cluster_num\Fig_3_barrier_clust_cluster_2_to_10.mat')
for j=1:length(foldername)
    load([foldername{j},'\','neuronIndividuals_new.mat'])
    for j1=7 % brick's empty
        for k=2:10
            [~,~,avg_region{j,j1}(k),~,~,~,avg_region_shuffled{j,j1}{k}]=DBSCAN_region_quantify_022422(group_ori_AI163_k{j,j1}{k},neuronIndividuals_new,[]);
        end
    end
end

% illustrate
for i=1:5
    for j=2:length(avg_region_shuffled{i,7})
        avg_region_shuffled{i,7}{j}=mean(avg_region_shuffled{i,7}{j});
    end
    avg_region_shuffled{i,7}=cell2mat(avg_region_shuffled{i,7});
end

figure;
for i=1:5
    subplot(1,5,i);
    plot(avg_region{i,7}(2:end));
    hold on;
    plot(avg_region_shuffled{i,7});

    disp(signrank(avg_region{i,7}(2:end),avg_region_shuffled{i,7}));
end

% statistical test
