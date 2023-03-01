function [group_period,A_color,A_color_region]=cluster_stability_illustrate_func(time_series,session,neuronIndividuals_new,win_leng,group_in)

total_leng=size(neuronIndividuals_new{session}.C,2);

%% 1. cluster num determine
if ~isempty(group_in)
    clust_size=length(unique(group_in));
    group_in_t{1}=group_in;
else
    disp('determining overall cluster size');
    [group_in_t]=cluster_determine_by_suoqin_NMF_firstPeakCoph(neuronIndividuals_new{session},100,10,[]);
    clust_size=length(unique(group_in_t));
end
clc
disp(['overall cluster size: ',num2str(clust_size)]);
group_in=group_in_t{1};

%% 2. calculate clusters across periods
ctt=1;
group_period={};

for i=time_series
    neuronIndividuals_temp=neuronIndividuals_new{session}.copy;
    neuronIndividuals_temp.C=neuronIndividuals_temp.C(:,i:min(i+win_leng-1,total_leng));
    neuronIndividuals_temp.S=neuronIndividuals_temp.S(:,i:min(i+win_leng-1,total_leng));
    
    [group,CM,Z,cophList,cophList_adjusted]=cluster_determine_by_suoqin_NMF_firstPeakCoph_stability(neuronIndividuals_temp,100,10,clust_size);
    group_period{ctt}=group;
    

    ctt=ctt+1;
end

%% draw anatomical map
%find corresponding group in current cluster result and assign its
A_color={};
A_color_region={};

for i=1:length(group_period)
    if ~isempty(group_in)
        % reindex calculated clusters so they looks as similar s possible
        [~,group_reindexed,~]=determineSharedCells_new(group_in,group_period{i}); 
        group_period{i}=group_reindexed;
    else
        if i>1
            % reindex calculated clusters so they looks as similar s possible
            [~,group_reindexed,~]=determineSharedCells_new(group_period{i-1},group_period{i});
            group_period{i}=group_reindexed;
        end
    end
    [A_color{i},A_color_region{i}]=DBSCAN_region_quantify_022422(group_period{i},neuronIndividuals_new,[]);
end

