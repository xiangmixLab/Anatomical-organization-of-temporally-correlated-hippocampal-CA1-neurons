% cross_time_stability_arranged
function [gp_all,A_color,A_color_region]=cross_time_stability_clust_cal(foldername,gp_trial)

% 1. calculate cluster
gp_all={};
neuron_all={};
for tk=1:length(foldername)
    load([foldername{tk},'\','neuronIndividuals_new.mat']);
    for j=1:length(neuronIndividuals_new)
        max_clust=max(gp_trial{tk,j});
        neuronCurr=neuronIndividuals_new{j}.copy;
        neuronIndividuals_cross_time=neuronIndividuals_new_split(neuronIndividuals_new{j},15*60*6);
        for k=1:length(neuronIndividuals_cross_time)
            [~,gp_all{tk,j}{k}]=cluster_determine_by_suoqin_NMF_firstPeakCoph_022422(neuronIndividuals_cross_time{k},100,10,max_clust);
            neuron_all{tk,j}{k}=neuronIndividuals_cross_time{k}.copy;
        end
    end
end

% 2. align clusters
for tk=1:size(gp_all,1)
    for j=1:size(gp_all,2)
        for k=2:length(gp_all{tk,j})
            [~,gp_all{tk,j}{k}]=determineSharedCells_new(gp_all{tk,j}{1},gp_all{tk,j}{k});
        end
    end
end

% 3. colored footprints and regions
A_color={};
A_color_region={};
for tk=1:size(gp_all,1)
    for j=1:size(gp_all,2)
        for k=1:length(gp_all{tk,j})
            [A_color{tk,j}{k},A_color_region{tk,j}{k}]=DBSCAN_region_quantify_022422(gp_all{tk,j}{k},{neuron_all{tk,j}{k}},[]);
        end
    end
end
