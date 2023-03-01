%% supplemental Fig. 8: Ai163 region analysis

foldername={
    'F:\052721_CA1 Vector Trace_Qiao_2021\101F\results'
    'F:\052721_CA1 Vector Trace_Qiao_2021\102F\results'
    'F:\052721_CA1 Vector Trace_Qiao_2021\103F\results'
    'F:\052721_CA1 Vector Trace_Qiao_2021\840F\results'
    'F:\052721_CA1 Vector Trace_Qiao_2021\2833M\results'
}


A_color={};
A_color_region={};
avg_region={};
max_region={};
avg_region_nnum={};
max_region_nnum={};
avg_region_shuf={};
max_region_shuf={};
avg_region_nnum_shuf={};
max_region_nnum_shuf={};
A_color_shuf={};
A_color_region_shuf={};
nd_all={};

tic;
for i=1:5
    load([foldername{i},'\','neuronIndividuals_new.mat']);
    for j=[5,8]
        for k=2:10
            group_AI163_regionNum{i,j}{k}=cluster_determine_by_suoqin_NMF_firstPeakCoph_022422(neuronIndividuals_new{j},100,10,k);
            [A_color{i,j}{k},A_color_region{i,j}{k},avg_region{i,j}{k},max_region{i,j}{k},avg_region_nnum{i,j}{k},max_region_nnum{i,j}{k},avg_region_shuf{i,j}{k},max_region_shuf{i,j}{k},avg_region_nnum_shuf{i,j}{k},max_region_nnum_shuf{i,j}{k},nd_all{i,j}(k),A_color_shuf{i,j}{k},A_color_region_shuf{i,j}{k}]=DBSCAN_region_quantify_022422(group_AI163_regionNum{i,j}{k},neuronIndividuals_new,[]);
            toc;
        end
    end
end

% 35um
tic;
for i=1:5
    load([foldername{i},'\','neuronIndividuals_new.mat']);
    for j=1:7
        for k=2:10
            [group_AI163_regionNum{i,j}{k}]=cluster_determine_by_suoqin_NMF_firstPeakCoph_022422(neuronIndividuals_new{j},100,10,k);
            [A_color{i,j}{k},A_color_region{i,j}{k},avg_region{i,j}{k},max_region{i,j}{k},avg_region_nnum{i,j}{k},max_region_nnum{i,j}{k},avg_region_shuf{i,j}{k},max_region_shuf{i,j}{k},avg_region_nnum_shuf{i,j}{k},max_region_nnum_shuf{i,j}{k},nd_all{i,j}(k),A_color_shuf{i,j}{k},A_color_region_shuf{i,j}{k}]=DBSCAN_region_quantify_022422(group_AI163_regionNum{i,j}{k},neuronIndividuals_new,17.5);
            toc;
        end
    end
end

%% illustrate behavior trace
load('F:\052721_CA1 Vector Trace_Qiao_2021\103F\results\training_cue_Behav.mat')
plot(behav.position(:,1),behav.position(:,2),'lineWidth',2);

%% intra-inter cluster corr
load('D:\Xu_clusterting_paper_prep11_2020\Round20\figS13_panels\AI163_clustering_ICA_kmean_originalExtract_040222_cond1-9.mat')
for tk=1:length(foldername)
    load([foldername{tk},'\neuronIndividuals_new.mat']);
    [intra_all,inter_all,~,~,intra_shuffle_all,~,~,intra_all_neuron{tk},inter_all_neuron{tk},intra_shuffle_all_neuron{tk}]=intra_inter_cluster_corr_dis(neuronIndividuals_new,group_ori_AI163{tk,9},9,'corr');
end

colorss=colormap(lines);
colorss=[colorss(4:6,:)];
cluster_correlation_cdf_plot({{intra_all},{inter_all},{intra_shuffle_all}},colorss);
[v,p]=intra_inter_statistic_perNeuron({intra_all_neuron,inter_all_neuron,intra_shuffle_all_neuron},{'intra','inter','intra_shuffle'});


%% calcualte region 
load('D:\Xu_clusterting_paper_prep11_2020\Round19\figS12_panels\AI163_clusterSize_data_95percentileRadius_originalExtract_033022_cond1_7.mat')
cond=9;

% current data is calculated as removing all neurons with max corr<=0.3
avg_region_cond1_mat=[];
for i=1:size(avg_region,1)
    avg_region_cond1_mat(i,:)=cell2mat(avg_region{i,cond}); 
end

avg_region_cond1_shuf_mat=[];
for i=1:size(avg_region_shuf,1) 
    shuf_region=avg_region_shuf{i,cond};
    shuf_region_mat=[];
    for j=1:length(shuf_region)
        shuf_region_mat(j)=mean(shuf_region{j});
    end
    avg_region_cond1_shuf_mat(i,:)=shuf_region_mat;
end

avg_region_cond1_shuf_90_mat=[];
for i=1:size(avg_region_shuf,1) 
    shuf_region=avg_region_shuf{i,cond};
    shuf_region_mat=[];
    for j=1:length(shuf_region)
        shuf_region_mat(j)=quantile(shuf_region{j},0.95);
    end
    avg_region_cond1_shuf_percentile_mat(i,:)=shuf_region_mat;
end

% illustrate single mice
for i=1:5
    subplot(1,5,i)
    plot(avg_region_cond1_mat(i,:));hold on;
    plot(avg_region_cond1_shuf_percentile_mat(i,2:end));
    plot(avg_region_cond1_shuf_mat(i,2:end));
    signrank(avg_region_cond1_mat(i,:),avg_region_cond1_shuf_mat(i,2:end))
end

% illustrate all mice with error bar, ALL
errorbar(mean(avg_region_cond1_mat(:,:),1),std(avg_region_cond1_mat(:,:),[],1)/size(avg_region_cond1_mat,1).^0.5)
hold on;
errorbar(mean(avg_region_cond1_shuf_percentile_mat(:,2:end),1),std(avg_region_cond1_shuf_percentile_mat(:,2:end),[],1)/(size(avg_region_cond1_shuf_percentile_mat(:,2:end),1)^0.5))
errorbar(mean(avg_region_cond1_shuf_mat(:,2:end),1),std(avg_region_cond1_shuf_mat(:,2:end)/size(avg_region_cond1_mat,1).^0.5,[],1))


% illustrate all mice with error bar, without 840F
errorbar(mean(avg_region_cond1_mat(1:4,:),1),std(avg_region_cond1_mat(1:4,:),[],1))
hold on;
errorbar(mean(avg_region_cond1_shuf_mat(1:4,2:end),1),std(avg_region_cond1_shuf_mat(1:4,2:end),[],1))

% current data is calculated as removing all neurons with max corr<=0.3
max_region_cond1_mat=[];
for i=1:size(max_region,1)
    max_region_cond1_mat(i,:)=cell2mat(max_region{i,cond});
end

max_region_cond1_shuf_mat=[];
for i=1:size(max_region_shuf,1)
    shuf_region=max_region_shuf{i,cond};
    shuf_region_mat=[];
    for j=1:length(shuf_region)
        shuf_region_mat(j)=mean(shuf_region{j});
    end
    max_region_cond1_shuf_mat(i,:)=shuf_region_mat;
end

errorbar(mean(max_region_cond1_mat,1),std(max_region_cond1_mat,[],1));
hold on;
errorbar(mean(max_region_cond1_shuf_mat(:,2:end),1),std(max_region_cond1_shuf_mat(:,2:end),[],1));

% current data is calculated as removing all neurons with max corr<=0.3
avg_region_nnum_cond1_mat=[];
for i=1:size(avg_region_nnum,1)
    avg_region_nnum_cond1_mat(i,:)=ceil(cell2mat(avg_region_nnum{i,1})); 
end

avg_region_nnum_cond1_shuf_mat=[];
for i=1:size(avg_region_nnum_shuf,1) 
    shuf_region_nnum=avg_region_nnum_shuf{i,1};
    shuf_region_nnum_mat=[];
    for j=1:length(shuf_region)
        shuf_region_nnum_mat(j)=mean(shuf_region_nnum{j});
    end
    avg_region_nnum_cond1_shuf_mat(i,:)=ceil(shuf_region_nnum_mat);
end

% 840F neuron density is smaller than the other four
errorbar(mean(avg_region_nnum_cond1_mat(:,:),1),std(avg_region_nnum_cond1_mat(:,:),[],1)/(size(avg_region_nnum_cond1_mat(:,:),1)^0.5))
hold on;
errorbar(mean(avg_region_nnum_cond1_shuf_mat(:,2:end),1),std(avg_region_nnum_cond1_shuf_mat(:,2:end),[],1)/(size(avg_region_nnum_cond1_shuf_mat(:,2:end),1)^0.5))

%% illustration
cond=6
subplot(121)
imagesc(A_color{3,cond}{3})
subplot(122)
imagesc(A_color_region{3,cond}{3})

figure;
ctt=1;
for i=1:10:100
    subplot(2,5,ctt);
    imagesc(A_color_shuf{3,cond}{3}{i});
    ctt=ctt+1;
end

figure;
ctt=1;
for i=1:10:100
    subplot(2,5,ctt);
    imagesc(A_color_region_shuf{3,cond}{3}{i});
    ctt=ctt+1;
end

%% compare size with optimal cluster number for each animal
load('D:\Xu_clusterting_paper_prep11_2020\Round19\figS12_panels\AI163_clustering_ICA_kmean.mat');
region_compare=[];
for i=1:size(group_ori_AI163,1)
    region_compare(i,:)=[avg_region_cond1_mat(i,max(group_ori_AI163{i,1})-1),avg_region_cond1_shuf_mat(i,max(group_ori_AI163{i,1}))];
end

max_region_compare=[];
for i=1:size(group_ori_AI163,1)
    max_region_compare(i,:)=[max_region_cond1_mat(i,max(group_ori_AI163{i,1})-1),max_region_cond1_shuf_mat(i,max(group_ori_AI163{i,1}))];
end

%% validation: multiGeo

foldername_multiGeo={
    'D:\Remapping_square_circle_triangle_061119_061319\M3411'	
    'D:\Remapping_square_circle_triangle_061119_061319\M3412'	
    'D:\Remapping_square_circle_triangle_061119_061319\M3421F'	
    'D:\Remapping_square_circle_triangle_061119_061319\M3422F'	
    'D:\Remapping_square_circle_triangle_061119_061319\M3424F'	
    'D:\Remapping_square_circle_triangle_061119_061319\M3425F'	
    }



% load('AI163_clustering_ICA_kmean.mat');

A_color_mGeo={};
A_color_region_mGeo={};
avg_region_mGeo={};
max_region_mGeo={};
avg_region_nnum_mGeo={};
max_region_nnum_mGeo={};
avg_region_shuf_mGeo={};
max_region_shuf_mGeo={};
avg_region_nnum_shuf_mGeo={};
max_region_nnum_shuf_mGeo={};
A_color_shuf_mGeo={};
A_color_region_shuf_mGeo={};
nd_all_mGeo={};

tic;
for i=1:6
    load([foldername_multiGeo{i},'\','neuronIndividuals_new.mat']);
    for j=3
        for k=2:10
            group_mGeo_regionNum{i,j}{k}=cluster_determine_by_suoqin_NMF_firstPeakCoph_022422(neuronIndividuals_new{j},100,10,k);
            [A_color_mGeo{i,j}{k},A_color_region_mGeo{i,j}{k},avg_region_mGeo{i,j}{k},max_region_mGeo{i,j}{k},avg_region_nnum_mGeo{i,j}{k},max_region_nnum_mGeo{i,j}{k},avg_region_shuf_mGeo{i,j}{k},max_region_shuf_mGeo{i,j}{k},avg_region_nnum_shuf_mGeo{i,j}{k},max_region_nnum_shuf_mGeo{i,j}{k},nd_all_mGeo{i,j}(k),A_color_shuf_mGeo{i,j}{k},A_color_region_shuf_mGeo{i,j}{k}]=DBSCAN_region_quantify_022422(group_mGeo_regionNum{i,j}{k},neuronIndividuals_new,[]);
            toc;
        end
    end
end

%%%%%%%%%%%
avg_region_cond1_mGeo_mat=[];
for i=1:size(avg_region_mGeo,1)
    avg_region_cond1_mGeo_mat(i,:)=cell2mat(avg_region_mGeo{i,3}); 
end

avg_region_cond1_shuf_mat=[];
for i=1:size(avg_region_shuf_mGeo,1) 
    shuf_region=avg_region_shuf_mGeo{i,3};
    shuf_region_mat=[];
    for j=1:length(shuf_region)
        shuf_region_mat(j)=mean(shuf_region{j});
    end
    avg_region_cond1_mGeo_shuf_mat(i,:)=shuf_region_mat;
end

%  neuron density is smaller than the other four
errorbar(mean(avg_region_cond1_mGeo_mat(1:4,:),1),std(avg_region_cond1_mGeo_mat(1:4,:),[],1)/(size(avg_region_cond1_mGeo_mat(1:4,:),1)^0.5))
hold on;
errorbar(mean(avg_region_cond1_mGeo_shuf_mat(1:4,2:end),1),std(avg_region_cond1_mGeo_shuf_mat(1:4,2:end),[],1)/(size(avg_region_cond1_mGeo_shuf_mat(1:4,2:end),1)^0.5))

% current data is calculated as removing all neurons with max corr<=0.3
max_region_cond1_mGeo_mat=[];
for i=1:size(max_region_mGeo,1)
    max_region_cond1_mGeo_mat(i,:)=cell2mat(max_region_mGeo{i,3});
end

max_region_cond1_mGeo_shuf_mat=[];
for i=1:size(max_region_shuf_mGeo,1)
    shuf_region=max_region_shuf_mGeo{i,3};
    shuf_region_mat=[];
    for j=1:length(shuf_region)
        shuf_region_mat(j)=mean(shuf_region{j});
    end
    max_region_cond1_mGeo_shuf_mat(i,:)=shuf_region_mat;
end

errorbar(mean(max_region_cond1_mGeo_mat,1),std(max_region_cond1_mGeo_mat,[],1)/(size(max_region_cond1_mGeo_mat,1)^0.5))
hold on;
errorbar(mean(max_region_cond1_mGeo_shuf_mat(:,2:end),1),std(max_region_cond1_mGeo_shuf_mat(:,2:end),[],1)/(size(max_region_cond1_mGeo_shuf_mat(:,2:end),1)^0.5))

% current data is calculated as removing all neurons with max corr<=0.3
avg_region_nnum_cond1_mGeo_mat=[];
for i=1:size(avg_region_nnum_mGeo,1)
    avg_region_nnum_cond1_mGeo_mat(i,:)=ceil(cell2mat(avg_region_nnum_mGeo{i,3})); 
end

avg_region_nnum_cond1_mGeo_shuf_mat=[];
for i=1:size(avg_region_nnum_shuf_mGeo,1) 
    shuf_region_nnum=avg_region_nnum_shuf_mGeo{i,3};
    shuf_region_nnum_mat=[];
    for j=1:length(shuf_region)
        shuf_region_nnum_mat(j)=mean(shuf_region_nnum{j});
    end
    avg_region_nnum_cond1_mGeo_shuf_mat(i,:)=ceil(shuf_region_nnum_mat);
end

%  neuron density is smaller than the other four
errorbar(mean(avg_region_nnum_cond1_mGeo_mat(1:4,:),1),std(avg_region_nnum_cond1_mGeo_mat(1:4,:),[],1)/(size(avg_region_nnum_cond1_mGeo_mat(1:4,:),1)^0.5))
hold on;
errorbar(mean(avg_region_nnum_cond1_mGeo_shuf_mat(1:4,2:end),1),std(avg_region_nnum_cond1_mGeo_shuf_mat(1:4,2:end),[],1)/(size(avg_region_nnum_cond1_mGeo_shuf_mat(1:4,2:end),1)^0.5))

% multiGeo illustration
subplot(121)
imagesc(A_color_mGeo{3,3}{5})
subplot(122)
imagesc(A_color_region_mGeo{3,3}{5})

figure;
ctt=1;
for i=1:10:100
    subplot(2,5,ctt);
    imagesc(A_color_shuf_mGeo{3,3}{5}{i});
    ctt=ctt+1;
end

figure;
ctt=1;
for i=1:10:100
    subplot(2,5,ctt);
    imagesc(A_color_region_shuf_mGeo{3,3}{5}{i});
    ctt=ctt+1;
end

%% original extraction clustering analysis
foldername={
    'F:\052721_CA1 Vector Trace_Qiao_2021\101F\results'
    'F:\052721_CA1 Vector Trace_Qiao_2021\102F\results'
    'F:\052721_CA1 Vector Trace_Qiao_2021\103F\results'
    'F:\052721_CA1 Vector Trace_Qiao_2021\840F\results'
    'F:\052721_CA1 Vector Trace_Qiao_2021\2833F\results'
}


A_color={};
A_color_region={};
avg_region={};
max_region={};
avg_region_nnum={};
max_region_nnum={};
avg_region_shuf={};
max_region_shuf={};
avg_region_nnum_shuf={};
max_region_nnum_shuf={};
A_color_shuf={};
A_color_region_shuf={};
nd_all={};

tic;
for i=1:5
    load([foldername{i},'\','neuronIndividuals_new.mat']);
    for j=[1,2,3,4,6,7,9]
        for k=2:10
            group_AI163_ori_regionNum{i,j}{k}=cluster_determine_by_suoqin_NMF_firstPeakCoph_022422(neuronIndividuals_new{j},100,10,k);
            [A_color{i,j}{k},A_color_region{i,j}{k},avg_region{i,j}{k},max_region{i,j}{k},avg_region_nnum{i,j}{k},max_region_nnum{i,j}{k},avg_region_shuf{i,j}{k},max_region_shuf{i,j}{k},avg_region_nnum_shuf{i,j}{k},max_region_nnum_shuf{i,j}{k},nd_all{i,j}(k),A_color_shuf{i,j}{k},A_color_region_shuf{i,j}{k}]=DBSCAN_region_quantify_022422(group_AI163_ori_regionNum{i,j}{k},neuronIndividuals_new,[]);
            toc;
        end
    end
end

%% original extraction, optimal clusters
load('D:\Xu_clusterting_paper_prep11_2020\Round19\figS13_panels\AI163_clustering_ICA_kmean_originalExtract_040222_cond1-9.mat');
for i=1:5
    foldername_AI163{i}(1:3)='F:\';
end
tic;
for i=1:5
    load([foldername_AI163{i},'\','neuronIndividuals_new.mat']);
    for j=9
        [~,~,avg_region{i,j},max_region{i,j},avg_region_nnum{i,j},max_region_nnum{i,j},avg_region_shuf{i,j},max_region_shuf{i,j},avg_region_nnum_shuf{i,j},max_region_nnum_shuf{i,j},nd_all{i,j},A_color_shuf{i,j},A_color_region_shuf{i,j}]=DBSCAN_region_quantify_022422(group_ori_AI163{i,j},neuronIndividuals_new,[]);
        toc;
    end
end
cond=9;
avg_region_mat=cell2mat(avg_region(:,cond));
avg_region_shuf_mat=[];
for i=1:size(avg_region_shuf,1) 
    shuf_region=avg_region_shuf{i,cond};
    avg_region_shuf_mat(i,1)=mean(shuf_region);
end
%% virus experiment 1, optimal clusters
load('D:\Xu_clusterting_paper_prep11_2020\Round19\figS13_panels\multiGeo_clustering_ICA_kmean.mat');
cond=3;

foldername_multiGeo={
    'D:\Remapping_square_circle_triangle_061119_061319\M3411'	
    'D:\Remapping_square_circle_triangle_061119_061319\M3412'	
    'D:\Remapping_square_circle_triangle_061119_061319\M3421F'	
    'D:\Remapping_square_circle_triangle_061119_061319\M3422F'	
    'D:\Remapping_square_circle_triangle_061119_061319\M3424F'	
    'D:\Remapping_square_circle_triangle_061119_061319\M3425F'	
    }

tic;
for i=1:6
    load([foldername_multiGeo{i},'\','neuronIndividuals_new.mat']);
    for j=1:6
        [~,~,avg_region{i,j},max_region{i,j},avg_region_nnum{i,j},max_region_nnum{i,j},avg_region_shuf{i,j},max_region_shuf{i,j},avg_region_nnum_shuf{i,j},max_region_nnum_shuf{i,j},nd_all{i,j},A_color_shuf{i,j},A_color_region_shuf{i,j}]=DBSCAN_region_quantify_022422(group_ori_multiGeo{i,j},neuronIndividuals_new,[]);
        toc;
    end
end

avg_region_mat=cell2mat(avg_region);
avg_region_shuf_mat=[];
for i=1:size(avg_region_shuf,1) 
    for cond=1:size(avg_region_shuf,2) 
        avg_region_shuf_mat(i,cond)=mean(avg_region_shuf{i,cond});
    end
end
