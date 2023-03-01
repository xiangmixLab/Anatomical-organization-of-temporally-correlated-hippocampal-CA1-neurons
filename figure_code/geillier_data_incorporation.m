% suppl 23. Geiller et al. data incorporation

% fig.4d: distance to seed neurons
% 1. get all pc: linear track, multiGeo, fig2, ai163
color_clust=distinguishable_colors(20);
load('D:\Xu_clusterting_paper_prep11_2020\Round22\figS20_panels\pc.mat');

% 2. get all clust
load('D:\Xu_clusterting_paper_prep11_2020\Round22\final_cluster_Data\multiGeo_clust.mat');
load('D:\Xu_clusterting_paper_prep11_2020\Round22\final_cluster_Data\fig2_clust.mat');
load('D:\Xu_clusterting_paper_prep11_2020\Round22\final_cluster_Data\AI163_clust.mat');
load('D:\Xu_clusterting_paper_prep11_2020\Round22\final_cluster_Data\LT_group.mat');

dis2clustCenThresh=1000; % some pcs are at edge, which may not those losonsky found, but still included
binw=20
% multiGeo
figure;
for j=1:size(all_pc_multiGeo,2)

    dis_bin_all={};
    pc_idx=[];
    count=1;
    
    for i=1:length(foldername_multiGeo)
        load([foldername_multiGeo{i},'\','neuronIndividuals_new.mat']);
        currPcs=all_pc_multiGeo{i,j}{2};
        currCen=neuronIndividuals_new{j}.centroid;
        groupAnatomical=major_anatomical_clusters(group_ori_multiGeo{i,j},currCen);
        for k=1:length(currPcs)
            if groupAnatomical(currPcs(k))>0
                currCenPC=currCen(currPcs(k),:);
                currCenClust=currCen(groupAnatomical==groupAnatomical(currPcs(k)),:);
                currCenClustcen=mean(currCenClust,1);
                
                if sum((currCenClustcen-currCenPC).^2,2)^0.5<=dis2clustCenThresh
                    dis=sum((currCenClust-currCenPC).^2,2).^0.5;
                    dis(dis==0)=[];
                    dis_bin=histcounts(dis*2,[0:binw:310]);

                    dis_bin_all{count,1}=dis_bin;
                    count=count+1;
                    pc_idx=[pc_idx;currPcs(k)];
                end
            end
        end
    end
    
    subplot(2,size(all_pc_multiGeo,2)/2,j);
    dis_bin_all_mat=cell2mat(dis_bin_all);
    dis_bin_all_mean=mean(dis_bin_all_mat,1);
    dis_bin_all_sem=sem(dis_bin_all_mat,1);

    errorbar([0:binw:binw*(length(dis_bin_all_mean)-1)],dis_bin_all_mean,dis_bin_all_sem)
    ylim([0,round(max(dis_bin_all_mean))+1])

end


% methology illustration
[d1,d2]=size(neuronIndividuals_new{1}.Cn);
subplot(131)
A_color=cluster_spatial_footprint_colormap(neuronIndividuals_new,d1,d2,color_clust,group_ori_multiGeo{i,j},0.75);
imagesc(A_color)

subplot(132)
groupAnatomical(groupAnatomical<=0)=-1;
A_color=cluster_spatial_footprint_colormap(neuronIndividuals_new,d1,d2,color_clust,groupAnatomical,0.75);
imagesc(A_color)

subplot(133)
gp=groupAnatomical*0;
gp(pc_idx)=1;
gp(gp==0)=-1;
A_color=cluster_spatial_footprint_colormap(neuronIndividuals_new,d1,d2,color_clust,gp,0.75);
imagesc(A_color)

figure;
k=137
% 
% for k=1:277 % one example pc 
%    gpt=ones(length(groupAnatomical),1)*-1;
%     gpt(group_ori_multiGeo{i,j}==2)=2;
%     gpt(k)=3;
%     A_color=cluster_spatial_footprint_colormap(neuronIndividuals_new,d1,d2,color_clust,gpt,0.75);
%     imagesc(A_color)
%     title(num2str(k));
%     drawnow;
%     pause(0.1);
%     clf;
% end

currCenPC=currCen(k,:);
currCenClust=currCen(groupAnatomical==groupAnatomical(k),:);
dis=sum((currCenClust-currCenPC).^2,2).^0.5*2;
dis=round(dis/(max(dis)/5))+1;
gpt=ones(length(groupAnatomical),1)*-1;
gpt(groupAnatomical==groupAnatomical(k))=dis;
A_color=cluster_spatial_footprint_colormap(neuronIndividuals_new,d1,d2,color_clust,gpt,0.75);
imagesc(A_color)

% fig2
figure;
for j=1:size(all_pc_fig2,2)

    dis_bin_all={};
    pc_idx=[];
    count=1;
    
    for i=1:length(foldername_fig2)
        load([foldername_fig2{i},'\','neuronIndividuals_new.mat']);
        currPcs=all_pc_fig2{i,j}{2};
        currCen=neuronIndividuals_new{j}.centroid;
        groupAnatomical=major_anatomical_clusters(group_ori_fig2{i,j},currCen);
        for k=1:length(currPcs)
            if groupAnatomical(currPcs(k))>0
                currCenPC=currCen(currPcs(k),:);
                currCenClust=currCen(groupAnatomical==groupAnatomical(currPcs(k)),:);
                currCenClustcen=mean(currCenClust,1);
                
                if sum((currCenClustcen-currCenPC).^2,2)^0.5<=dis2clustCenThresh
                    dis=sum((currCenClust-currCenPC).^2,2).^0.5;
                    dis(dis==0)=[];
                    dis_bin=histcounts(dis*2,[0:binw:310]);

                    dis_bin_all{count,1}=dis_bin;
                    count=count+1;
                    pc_idx=[pc_idx;currPcs(k)];
                end
            end
        end
    end
    
    subplot(1,size(all_pc_fig2,2),j);
    dis_bin_all_mat=cell2mat(dis_bin_all);
    dis_bin_all_mean=mean(dis_bin_all_mat,1);
    dis_bin_all_sem=sem(dis_bin_all_mat,1);

    errorbar([0:binw:binw*(length(dis_bin_all_mean)-1)],dis_bin_all_mean,dis_bin_all_sem)
    ylim([0,round(max(dis_bin_all_mean))+1])
end

% blocker
figure;
for j=1:size(all_pc_AI163,2)

    dis_bin_all={};
    pc_idx=[];
    count=1;
    
    for i=1:length(foldername_AI163)
        load([foldername_AI163{i},'\','neuronIndividuals_new.mat']);
        currPcs=all_pc_AI163{i,j}{2};
        currCen=neuronIndividuals_new{j}.centroid;
        groupAnatomical=major_anatomical_clusters(group_ori_AI163{i,j},currCen);
        for k=1:length(currPcs)
            if groupAnatomical(currPcs(k))>0
                currCenPC=currCen(currPcs(k),:);
                currCenClust=currCen(groupAnatomical==groupAnatomical(currPcs(k)),:);
                currCenClustcen=mean(currCenClust,1);
                
                if sum((currCenClustcen-currCenPC).^2,2)^0.5<=dis2clustCenThresh
                    dis=sum((currCenClust-currCenPC).^2,2).^0.5;
                    dis(dis==0)=[];
                    dis_bin=histcounts(dis*2,[0:binw:310]);

                    dis_bin_all{count,1}=dis_bin;
                    count=count+1;
                    pc_idx=[pc_idx;currPcs(k)];
                end
            end
        end
    end
    
    subplot(2,size(all_pc_AI163,2)/2,j);
    dis_bin_all_mat=cell2mat(dis_bin_all);
    dis_bin_all_mean=mean(dis_bin_all_mat,1);
    dis_bin_all_sem=sem(dis_bin_all_mat,1);

    errorbar([0:binw:binw*(length(dis_bin_all_mean)-1)],dis_bin_all_mean,dis_bin_all_sem)
    ylim([0,round(max(dis_bin_all_mean))+1])
end

% linear track
figure;
for j=1:size(all_pc_lt,2)

    dis_bin_all={};
    pc_idx=[];
    count=1;
    
    for i=1:length(foldername_lt)
        load([foldername_lt{i},'\','neuronIndividuals_new.mat']);
        currPcs=all_pc_lt{i,j}{2};
        currCen=neuronIndividuals_new{j}.centroid;
        groupAnatomical=major_anatomical_clusters(group_ori_LT{i,j},currCen);
        for k=1:length(currPcs)
            if groupAnatomical(currPcs(k))>0
                currCenPC=currCen(currPcs(k),:);
                currCenClust=currCen(groupAnatomical==groupAnatomical(currPcs(k)),:);
                currCenClustcen=mean(currCenClust,1);
                
                if sum((currCenClustcen-currCenPC).^2,2)^0.5<=dis2clustCenThresh
                    dis=sum((currCenClust-currCenPC).^2,2).^0.5;
                    dis(dis==0)=[];
                    dis_bin=histcounts(dis*2,[0:binw:310]);

                    dis_bin_all{count,1}=dis_bin;
                    count=count+1;
                    pc_idx=[pc_idx;currPcs(k)];
                end
            end
        end
    end
    
    subplot(1,size(all_pc_lt,2),j);
    dis_bin_all_mat=cell2mat(dis_bin_all);
    dis_bin_all_mean=mean(dis_bin_all_mat,1);
    dis_bin_all_sem=sem(dis_bin_all_mat,1);

    errorbar([0:binw:binw*(length(dis_bin_all_mean)-1)],dis_bin_all_mean,dis_bin_all_sem)
    ylim([0,round(max(dis_bin_all_mean))+1])
end


% behav utility
behavname={
    'triangle1_Behav.mat';
    'circle1_Behav.mat'; 
    'square1_Behav.mat';
    'circle2_Behav.mat';
    'square2_Behav.mat';
    'triangle2_Behav.mat';
    }

all_behav={};
for i=1:length(foldername_multiGeo)
    for j=1:length(behavname)
        load([foldername_multiGeo{i},'\',behavname{j}]);
        all_behav{i,j}=behav;
        all_behav{i,j}.VidObj=[];
    end
end
