%DBSCAN region properities for different num of clusters
function [A_color,A_color_region,avg_max_region,avg_max_reg_shuf_all,avg_max_reg_shuf_all_95,A_color_region_shuffle,all_mean_region,avg_max_region_nnum,all_nnums,all_max_regions,nd,all_nnums_shuffled,region_all,region_all_shuffled,avg_max_region_nnum_shuffled]=DBSCAN_region_quantify_func_single_group_no_plot(group_record,neuronIndividuals_new,nd_given,recordingMethod)

clust_idx=length(unique(group_record));
shuffles=100;
neuronIndividuals_new{1}.A=full(neuronIndividuals_new{1}.A);
[boundary_all,region_all,region_nnum_all,nd]=cluster_region_cal(group_record,1,neuronIndividuals_new,nd_given,recordingMethod{1});
for k=1:shuffles
    r1=randperm(length(group_record));
    r2=r1(randperm(length(group_record)));
    group_rand= group_record(r2);
    [boundary_all_shuffled{k},region_all_shuffled{k},region_nnum_all_shuffled{k}]=cluster_region_cal(group_rand,1,neuronIndividuals_new,nd_given,recordingMethod{1});      
end

% 

% load(['C:\Users\exx\Desktop\HDAC paper fig and method\SFN2019\SFN2019 fig\DBSCAN regions\triangle_square_circle\M3412\boundary_and_region.mat']);

d1=neuronIndividuals_new{1}.imageSize(1);
d2=neuronIndividuals_new{1}.imageSize(2);
avg_max_region_shuffled={};
avg_max_peri_shuffled={};
area_thresh_all_shuffled={};
[avg_max_region,avg_max_peri,~,~,area_thresh_all,all_mean_region,all_max_regions,avg_max_region_nnum,all_nnums]=DBSCAN_region_perimeter_quantify_func_052821({region_all},{region_nnum_all});
for k=1:shuffles
    [avg_max_region_shuffled{k},avg_max_peri_shuffled{k},~,~,area_thresh_all_shuffled{k},all_mean_region_shuffled{k},all_max_regions_shuffled{k},avg_max_region_nnum_shuffled{k},all_nnums_shuffled{k}]=DBSCAN_region_perimeter_quantify_func_052821({region_all_shuffled{k}},{region_nnum_all_shuffled{k}});
end

avg_max_reg_shuf_all_dat=zeros(1,shuffles);
for p=1:shuffles
    avg_max_reg_shuf_all_dat(:,p)=avg_max_region_shuffled{p};
end
avg_max_reg_shuf_all=mean(avg_max_reg_shuf_all_dat,2);
avg_max_reg_shuf_all_95=quantile(avg_max_reg_shuf_all_dat,.95,2);


avg_max_peri_shuf_all_dat=zeros(1,shuffles);
for p=1:shuffles
    avg_max_peri_shuf_all_dat(1,p)=avg_max_peri_shuffled{p};
end
avg_max_peri_shuf_all=mean(avg_max_peri_shuf_all_dat,2);
avg_max_peri_shuf_all_95=quantile(avg_max_peri_shuf_all_dat,.95,2);

%clust_idx
if isempty(clust_idx)
    clust_idx=1;
end

%inneed: footprints
colorClusters_all=distinguishable_colors(100);
d1=neuronIndividuals_new{1}.imageSize(1);
d2=neuronIndividuals_new{1}.imageSize(2);
group=group_record; % just use the 4 cluster case as an example

thres_val=0.7;
if isequal(recordingMethod{1},'2p')
    thres_val=0.15;
end
A_color=cluster_spatial_footprint_colormap(neuronIndividuals_new,d1,d2,colorClusters_all,group,thres_val);

%inneed: regions
A_color_region=ones(d1*d2,3)*255;
for i=1:size(region_all,1)
    for j=1:size(region_all,2)
        if ~isempty(region_all{i,j})
            region_t=region_all{i,j}>0;
            region_t=bwareaopen(region_t,round(area_thresh_all{1}(i)));
            ind=find(region_t>0);%the pixels at which a neuron shows up
            if min(group)~=-1
                A_color_region(ind,:)=repmat(colorClusters_all(i,:),length(ind),1);
            else
                if i==1
                    A_color_region(ind,:)=repmat([0.8,0.8,0.8],length(ind),1);    
                else
                    A_color_region(ind,:)=repmat(colorClusters_all(i-1,:),length(ind),1);
                end
            end
            region_all{i,j}=region_t;
        end
    end
end
A_color_region=reshape(A_color_region,d1,d2,3);

%inneed: shuffle regions
% A_color_region_shuffle={};
% for k=1:length(region_all_shuffled)
%     A_color_region_shuffle_t=ones(d1*d2,3)*255;
%     for i=1:size(region_all_shuffled{k},1)
%         for j=1:size(region_all_shuffled{k},2)
%             if ~isempty(region_all_shuffled{k}{i,j}) 
%                 region_t=region_all_shuffled{k}{i,j}>0;
%                 try
%                     region_t=bwareaopen(region_t,round(area_thresh_all_shuffled{k}{1}(i)));
%                 catch
%                     region_t=bwareaopen(region_t,0);
%                 end
%                 ind=find(region_t>0);%the pixels at which a neuron shows up
%                 if min(group)~=-1
%                     A_color_region(ind,:)=repmat(colorClusters_all(i,:),length(ind),1); 
%                 else
%                     if i==1
%                         A_color_region(ind,:)=repmat([0.8,0.8,0.8],length(ind),1);    
%                     else
%                         A_color_region(ind,:)=repmat(colorClusters_all(i-1,:),length(ind),1);
%                     end
%                 end
%                 region_all_shuffled{k}{i,j}=region_t;
%             end
%         end
%     end
%     A_color_region_shuffle_t=reshape(A_color_region_shuffle_t,d1,d2,3);
%     A_color_region_shuffle{k}=A_color_region_shuffle_t;
    A_color_region_shuffle{k}=[];
end


