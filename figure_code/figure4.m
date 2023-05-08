%% immobility
%% actually, immobility. as mice will be immobile for a short period and may shift position in middle, for persentation only the first immobility session is shown here

behav_vid={
    'D:\Xu_clusterting_paper_prep11_2020\arranged_final_data\Fig_4_immobile\3422_result_102919\mobile.tif','D:\Xu_clusterting_paper_prep11_2020\arranged_final_data\Fig_4_immobile\3422_result_102919\immobile.tif';
    'D:\Xu_clusterting_paper_prep11_2020\arranged_final_data\Fig_4_immobile\3424_result_102919\mobile.tif','D:\Xu_clusterting_paper_prep11_2020\arranged_final_data\Fig_4_immobile\3424_result_102919\immobile.tif';
    'D:\Xu_clusterting_paper_prep11_2020\arranged_final_data\Fig_4_immobile\3425_result_102919\mobile.tif','D:\Xu_clusterting_paper_prep11_2020\arranged_final_data\Fig_4_immobile\3425_result_102919\immobile.tif';  
    };

neuronDat={
    'D:\Xu_clusterting_paper_prep11_2020\arranged_final_data\Fig_4_immobile\3422_result_102919\neuronIndividuals_new.mat';
    'D:\Xu_clusterting_paper_prep11_2020\arranged_final_data\Fig_4_immobile\3424_result_102919\neuronIndividuals_new.mat';
    'D:\Xu_clusterting_paper_prep11_2020\arranged_final_data\Fig_4_immobile\3425_result_102919\neuronIndividuals_new.mat';
    }

behavDat={
    'D:\Xu_clusterting_paper_prep11_2020\arranged_final_data\Fig_4_immobile\3422_result_102919\behav.mat';
    'D:\Xu_clusterting_paper_prep11_2020\arranged_final_data\Fig_4_immobile\3424_result_102919\behav.mat';
    'D:\Xu_clusterting_paper_prep11_2020\arranged_final_data\Fig_4_immobile\3425_result_102919\behav.mat';
    }

clustDat='D:\Xu_clusterting_paper_prep11_2020\final_code\final_cluster_data\cluster_optimal_num\Fig4_imb_clust_original.mat';
clustDat_2_10='D:\Xu_clusterting_paper_prep11_2020\final_code\final_cluster_data\cluster_num_2_to_10\Fig4_imb_clust_cluster_2_to_10.mat';


load(clustDat);
load(clustDat_2_10);

for i=1:3
    figure;
    set(gcf,'renderer','painter');    
    
    group_1=[];
    for j=1:size(behav_vid,2)
        
        % load behav vid
        frame=imread(behav_vid{i,j});
        
        % load neuron
        load(neuronDat{i});
        ntemp=neuronIndividuals_new{j};
        thresh=3*std(ntemp.S,[],2);
        ntemp.S=thresh_nC_nS(ntemp.S,thresh);
        
        % load behav
        load(behavDat{i});
        behavpos=behavIndividuals{j}.position;
        behavtime=behavIndividuals{j}.time;
        ntemp.time=behavtime;
        btime=behavtime;
        behav=behavIndividuals{j};
        
        % pickup first 5 min recording (make the time the same for mobile
        % and dark-immobile)
        time_cut=min(min(9000,size(ntemp.C,2)),size(behavpos,1))-1; %frame, 30Hz, first 5 min

        d1=240;
        d2=376;
        ntemp.C=ntemp.C(:,1:time_cut/2);
        ntemp.S=ntemp.S(:,1:time_cut/2);
        ntemp.time=ntemp.time(1:time_cut);
        btime=btime(1:time_cut,:);
        behavpos=behavpos(1:time_cut,:);
                 
        
        % calculate footprint
        group=group_ori_imb{i,j};
        colorClusters_all=distinguishable_colors(10);
%         A_color=cluster_spatial_footprint_colormap({ntemp},d1,d2,colorClusters_all,group,0.7);
         
        % area size across different clust num, compared with baseline
        group_rec=group_ori_imb_k{i,j};

         tic;
         A_color={};
         A_color_region={};
         avg_region=[];
         avg_reg_shuf=[];
         avg_reg_shuf_95=[];
         for pk=2:10
            [A_color{pk},A_color_region{pk},avg_region(pk),~,~,~,avg_reg_shuf_temp]=DBSCAN_region_quantify_022422(group_rec{pk},{ntemp},[]);
            avg_reg_shuf(pk)=mean(avg_reg_shuf_temp);
         end
         toc;
         area_size_summary{i,j}={avg_region,avg_reg_shuf};
         
         clust_idx=length(unique(group))-1;      
         
         %frame and traj
         subplot(size(behav_vid,2),6,1+(j-1)*6);
         imagesc(frame);
         hold on;
         plot(behavpos(:,1)*behav.ROI3/behav.trackLength+behav.ROI(1),behavpos(:,2)*behav.ROI3/behav.trackLength+behav.ROI(2),'color','c');
         
         
         %calcium heatmap
         subplot(size(behav_vid,2),6,2+(j-1)*6);
         dataC2 = [];
         for p = 1:length(unique(group))
             dataC2 = [dataC2;ntemp.C(group == p,:)];
         end

         positionC2 = 0;
         for p =1:length(unique(group))
             positionC2(p+1) = positionC2(p)+sum(group == p);
         end
         optimalK2=length(unique(group));

         imagesc(dataC2)
         colormap(jet)
         colorbar
         caxis([0 100]);
         hold on
         for p = 2:length(positionC2)-1
             line(get(gca,'XLim'),[positionC2(p) positionC2(p)],'LineWidth',0.5,'Color','r'); hold on;
         end
         set(gca,'FontSize',8)
         labels = strcat('C',cellstr(num2str([1:optimalK2]')));
         ytick0 = positionC2(1:end-1)+diff(positionC2)/2;
         set(gca,'Ytick',ytick0);set(gca,'YtickLabel',labels,'FontName','Arial','FontSize',10)
         xlabel('Frames','FontSize',10)
         ylabel('Neurons','FontSize',10)

         group2=group;
         
         % pairwise correlation
         subplot(size(behav_vid,2),6,3+(j-1)*6);
         perm=[];
         for p = 1:length(unique(group2))
             perm = [perm;find(group2 == p)];
         end
         CM_in=squareform(1-pdist(ntemp.C,'correlation'));
%          displayHeatmap_general(CM_in,CM_in,perm,colorClusters2,optimalK2,group2,1)

         CM_out=clustPairwiseHeatmap_gen(CM_in,CM_in,perm,1);
         imagesc(CM_out);
         pcorr_summary{i,j}=CM_out;
         group_summary{i,j}=group2;
         colormap(jet);
         caxis([0 0.8]);
         colorbar
         %footprint
         subplot(size(behav_vid,2),6,4+(j-1)*6);
         imagesc(A_color{max(group2)});
         %regionclea
         subplot(size(behav_vid,2),6,5+(j-1)*6);
         imagesc(A_color_region{max(group2)});
         %region size change
         subplot(size(behav_vid,2),6,6+(j-1)*6);
         plot(avg_region(2:end),'-o','color',[0    0.4470    0.7410]);
         hold on;
         plot(avg_reg_shuf(2:end),'-s','color',[0.8500    0.3250    0.0980]);
         fr_temp=sum(ntemp.S>0,2)./size(ntemp.S,2);
         fr_temp=fr_temp(fr_temp>10/size(ntemp.S,2)); % remove some small firing rate neurons
         fr_summary{i,j}=mean(fr_temp);
         all_avg_region{i,j}=avg_region;
    end
end

fr_summary_mat=cell2mat(fr_summary);

%% plot intra/inter cluster temporal corr
subplot(121)

h1=cdfplot(cell2mat(intra_all(:,1)));
hold on;
h2=cdfplot(cell2mat(inter_all(:,1)));
h3=cdfplot(cell2mat(intra_shuffle_all(:,1)));
set(h1,'color',[0.4940,0.1840,0.5560]);
set(h2,'color',[0.4660,0.6740,0.1880]);
set(h3,'color',[0.3010,0.7450,0.9330]);

p1=infer_cdf_loc(cell2mat(intra_all(:,1)),mean(cell2mat(intra_all(:,1))));
p2=infer_cdf_loc(cell2mat(inter_all(:,1)),mean(cell2mat(inter_all(:,1))));
p3=infer_cdf_loc(cell2mat(intra_shuffle_all(:,1)),mean(cell2mat(intra_shuffle_all(:,1))));

plot(p1(1),p1(2),'.','color',colorClusters_all(4,:),'MarkerSize',20)
plot(p2(1),p2(2),'.','color',colorClusters_all(5,:),'MarkerSize',20)
plot(p3(1),p3(2),'.','color',colorClusters_all(6,:),'MarkerSize',20)

subplot(122)

h1=cdfplot(cell2mat(intra_all(:,2)));
hold on;
h2=cdfplot(cell2mat(inter_all(:,2)));
h3=cdfplot(cell2mat(intra_shuffle_all(:,2)));
set(h1,'color',[0.4940,0.1840,0.5560]);
set(h2,'color',[0.4660,0.6740,0.1880]);
set(h3,'color',[0.3010,0.7450,0.9330]);

p1=infer_cdf_loc(cell2mat(intra_all(:,2)),mean(cell2mat(intra_all(:,2))));
p2=infer_cdf_loc(cell2mat(inter_all(:,2)),mean(cell2mat(inter_all(:,2))));
p3=infer_cdf_loc(cell2mat(intra_shuffle_all(:,2)),mean(cell2mat(intra_shuffle_all(:,2))));

plot(p1(1),p1(2),'.','color',colorClusters_all(4,:),'MarkerSize',20)
plot(p2(1),p2(2),'.','color',colorClusters_all(5,:),'MarkerSize',20)
plot(p3(1),p3(2),'.','color',colorClusters_all(6,:),'MarkerSize',20)

mean(cell2mat(intra_all(:,1)))
std(cell2mat(intra_all(:,1)))/length(cell2mat(intra_all(:,1)))^0.5
mean(cell2mat(inter_all(:,1)))
std(cell2mat(inter_all(:,1)))/length(cell2mat(inter_all(:,1)))^0.5
mean(cell2mat(intra_shuffle_all(:,1)))
std(cell2mat(intra_shuffle_all(:,1)))/length(cell2mat(intra_shuffle_all(:,1)))^0.5

k1=cell2mat(intra_all(:,1));
k2=cell2mat(inter_all(:,1));
k3=cell2mat(intra_shuffle_all(:,1));

mean(cell2mat(intra_all(:,2)))
std(cell2mat(intra_all(:,2)))/length(cell2mat(intra_all(:,2)))^0.5
mean(cell2mat(inter_all(:,2)))
std(cell2mat(inter_all(:,2)))/length(cell2mat(inter_all(:,2)))^0.5
mean(cell2mat(intra_shuffle_all(:,2)))
std(cell2mat(intra_shuffle_all(:,2)))/length(cell2mat(intra_shuffle_all(:,2)))^0.5

k1=cell2mat(intra_all(:,2));
k2=cell2mat(inter_all(:,2));
k3=cell2mat(intra_shuffle_all(:,2));


