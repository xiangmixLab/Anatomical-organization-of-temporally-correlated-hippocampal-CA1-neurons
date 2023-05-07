%% sleep 
%% actually, immobility. as mice will be immobile for a short period and may shift position in middle, for persentation only the first immobility session is shown here
behav_vid={
    'X:\Lab\DataImages\Yanjun\Yanjun_nnRevision\Miniscope_OLM_raw\12_18_2018_Baseline\M3422F\behavCam2.avi','K:\Cluster_paper_data\Sleep exp\08082019_sleep_3422\3422_sleep1_7min_last10seconds wake\behavCam1_2.avi';
    'X:\Lab\DataImages\Yanjun\Yanjun_nnRevision\Miniscope_OLM_raw\12_18_2018_Baseline\M3424F\behavCam2.avi','K:\Cluster_paper_data\Sleep exp\08092019_sleep_3424\sleep1_8min_last1min wakeup\behavCam1_2.avi';
    'X:\Lab\DataImages\Yanjun\Yanjun_nnRevision\Miniscope_OLM_raw\12_18_2018_Baseline\M3425F\behavCam2.avi','K:\Cluster_paper_data\Sleep exp\08072019_Sleep_3425\sleep6_5min_not move_start move last10s\behavCam1_2.avi';  
    'X:\Lab\DataImages\Yanjun\Yanjun_nnRevision\Miniscope_OLM_raw\12_18_2018_Baseline\M3425F\behavCam2.avi','K:\Cluster_paper_data\Sleep exp\08072019_Sleep_3425\sleep6_5min_not move_start move last10s\behavCam1_2.avi';    
    'X:\Lab\DataImages\Yanjun\Yanjun_nnRevision\Miniscope_OLM_raw\12_18_2018_Baseline\M3425F\behavCam2.avi','K:\Cluster_paper_data\Sleep exp\08072019_Sleep_3425\sleep6_5min_not move_start move last10s\behavCam1_2.avi';    
    };


neuronDat={
    'K:\Cluster_paper_data\Sleep exp\3422_result_102919\sleep\further_processed_neuron_extraction_final_result.mat','K:\Cluster_paper_data\Sleep exp\3422_result_102919\sleep\further_processed_neuron_extraction_final_result.mat';
    'K:\Cluster_paper_data\Sleep exp\3424_result_102919\active_sleep\further_processed_neuron_extraction_final_result.mat','K:\Cluster_paper_data\Sleep exp\3424_result_102919\active_sleep\further_processed_neuron_extraction_final_result.mat';
    'K:\Cluster_paper_data\Sleep exp\3425_result_102919\further_processed_neuron_extraction_final_result.mat','K:\Cluster_paper_data\Sleep exp\3425_result_102919\further_processed_neuron_extraction_final_result.mat';
    'K:\Cluster_paper_data\Sleep exp\270_result\further_processed_neuron_extraction_final_result.mat','K:\Cluster_paper_data\Sleep exp\270_result\further_processed_neuron_extraction_final_result.mat';
    'K:\Cluster_paper_data\Sleep exp\2833_result\further_processed_neuron_extraction_final_result.mat','K:\Cluster_paper_data\Sleep exp\2833_result\further_processed_neuron_extraction_final_result.mat';
    }

behavDat={
    'K:\Cluster_paper_data\Sleep exp\3422_result_102919\sleep\awake_results\current_condition_behav.mat',{'K:\Cluster_paper_data\Sleep exp\3422_result_102919\sleep\sleep1_results\current_condition_behav.mat','K:\Cluster_paper_data\Sleep exp\3422_result_102919\sleep\sleep2_results\current_condition_behav.mat','K:\Cluster_paper_data\Sleep exp\3422_result_102919\sleep\sleep3_results\current_condition_behav.mat'};
    'K:\Cluster_paper_data\Sleep exp\3424_result_102919\active_sleep\active_results\current_condition_behav.mat',{'K:\Cluster_paper_data\Sleep exp\3424_result_102919\active_sleep\sleep1_results\current_condition_behav.mat','K:\Cluster_paper_data\Sleep exp\3424_result_102919\active_sleep\sleep2_results\current_condition_behav.mat','K:\Cluster_paper_data\Sleep exp\3424_result_102919\active_sleep\sleep3_results\current_condition_behav.mat'};
    'K:\Cluster_paper_data\Sleep exp\3425_result_102919\active_results\current_condition_behav.mat',{'K:\Cluster_paper_data\Sleep exp\3425_result_102919\sleep1_results\current_condition_behav.mat','K:\Cluster_paper_data\Sleep exp\3425_result_102919\sleep2_results\current_condition_behav.mat','K:\Cluster_paper_data\Sleep exp\3425_result_102919\sleep3_results\current_condition_behav.mat','K:\Cluster_paper_data\Sleep exp\3425_result_102919\sleep4_results\current_condition_behav.mat','K:\Cluster_paper_data\Sleep exp\3425_result_102919\sleep5_results\current_condition_behav.mat','K:\Cluster_paper_data\Sleep exp\3425_result_102919\sleep6_results\current_condition_behav.mat'};
    'K:\Cluster_paper_data\Sleep exp\270_result\active_bright_Behav.mat',{'K:\Cluster_paper_data\Sleep exp\270_result\sleep_dark_Behav.mat'};
    'K:\Cluster_paper_data\Sleep exp\2833_result\active_bright_Behav.mat',{'K:\Cluster_paper_data\Sleep exp\2833_result\sleep_dark_Behav.mat'};

    }

behavExp={
    'K:\Cluster_paper_data\Sleep exp\3422_result_102919\sleep\awake_Behav.mat','K:\Cluster_paper_data\Sleep exp\3422_result_102919\sleep\sleep1_Behav.mat';
    'K:\Cluster_paper_data\Sleep exp\3424_result_102919\active_sleep\active_Behav.mat','K:\Cluster_paper_data\Sleep exp\3424_result_102919\active_sleep\sleep1_Behav.mat';
    'K:\Cluster_paper_data\Sleep exp\3425_result_102919\active_Behav.mat','K:\Cluster_paper_data\Sleep exp\3425_result_102919\sleep6_Behav.mat';
    'K:\Cluster_paper_data\Sleep exp\270_result\active_bright_Behav.mat','K:\Cluster_paper_data\Sleep exp\270_result\sleep_dark_Behav.mat';
    'K:\Cluster_paper_data\Sleep exp\2833_result\active_bright_Behav.mat','K:\Cluster_paper_data\Sleep exp\2833_result\sleep_dark_Behav.mat';
    }


cond_idx=[
    1 2;
    1 2;
    1 7;
    1 2;
    1 2;];

original_area_size_summary={};
area_size_summary={};
pcorr_summary={};
group_summary={};
fr_summary={};
intra_all={};
inter_all={};
intra_shuffle_all={};
intra_all_dis={};
inter_all_dis={};
intra_shuffle_all_dis={};

velo_summary=[];
all_avg_region={};

for i=1:3
    figure;
    set(gcf,'renderer','painter');    
    group_1=[];
    for j=1:size(behav_vid,2)
        load(neuronDat{i,j});
        vid=VideoReader(behav_vid{i,j});
        frame=readFrame(vid);
     
        if j==1
            ntemp=neuron.copy;
            ntemp.C=ntemp.C(:,1:ntemp.num2read(2));
            ntemp.S=C_to_peakS(ntemp.C);
       
            thresh=3*std(ntemp.S,[],2);
            ntemp.S=thresh_nC_nS(ntemp.S,thresh);
            
            load(behavDat{i,j});
            dist_chg=[];
            for k=1:size(behavpos,1)-1
                dist_chg(k)=norm(behavpos(k+1,:)-behavpos(k,:));
            end
            ntemp.time=behavtime;
            btime=behavtime;
        else
            ntemp=neuron.copy;
            ntemp.C=ntemp.C(:,ntemp.num2read(2)+1:end);
            ntemp.S=C_to_peakS(ntemp.C);
            
            thresh=3*std(ntemp.S,[],2);
            ntemp.S=thresh_nC_nS(ntemp.S,thresh);
            
            behavpos_all=[];
            behavtime_all=[];
            for k=1:length(behavDat{i,j})
                load(behavDat{i,j}{k});
                behavpos_all=[behavpos_all;behavpos];
                behavtime_all=[behavtime_all;behavtime];
            end
            behavpos=behavpos_all;
            dist_chg=[];
            for k=1:size(behavpos,1)-1
                dist_chg(k)=norm(behavpos(k+1,:)-behavpos(k,:));
            end
            
            Fs=30;
%             moving_session=dist_chg>=1;
%             moving_session1=moving_session;
%             for k=2*round(Fs)+1:length(moving_session1)-2*round(Fs)
%                 if moving_session1(k)==1
%                     moving_session(k-2*round(Fs):k+2*round(Fs))=1; % adapted Hai's suggestion to make +-1s of single velocity point become a period
%                 end
%             end
%             
%             dist_chg=moving_session;
            time_resample_idx=dist_chg<1;
            c_resample_idx=resample(double(time_resample_idx),size(ntemp.C,2),size(behavpos,1));
            c_resample_idx=logical(c_resample_idx);
            ntemp.C=ntemp.C(:,c_resample_idx);
            ntemp.S=ntemp.S(:,c_resample_idx);
            ntemp.time=behavtime_all(time_resample_idx);
            btime=behavtime_all(time_resample_idx);
            behavpos=behavpos(time_resample_idx,:);

        end
        
        time_cut=min(min(9000,size(ntemp.C,2)),size(behavpos,1))-1; %frame, 30Hz, first 5 min

        dist_chg=dist_chg(1,1:time_cut);
        d1=ntemp.imageSize(1);
        d2=ntemp.imageSize(2);
        ntemp.C=ntemp.C(:,1:time_cut/2);
        ntemp.S=ntemp.S(:,1:time_cut/2);
        ntemp.time=ntemp.time(1:time_cut);
        btime=btime(1:time_cut,:);
        behavpos=behavpos(1:time_cut,:);

        if j==1 % force sleep and awake have the same number of clusters      
            [group,CMO,Z,cophList,cophList_adjusted]=cluster_determine_by_suoqin_NMF_firstPeakCoph_stability(ntemp,100,10,[]);
            [~,~,original_area_size_summary{i,j}{1},~,~,~,original_area_size_summary{i,j}{2},~]=DBSCAN_region_quantify_022422(group,{ntemp},[]); % make sleep clust # the same with awake
        else
            [group,CMO,Z,cophList,cophList_adjusted]=cluster_determine_by_suoqin_NMF_firstPeakCoph_stability(ntemp,100,10,max(group_summary{i,1}));
            [~,~,original_area_size_summary{i,j}{1},~,~,~,original_area_size_summary{i,j}{2},~]=DBSCAN_region_quantify_022422(group,{ntemp},[]); % make sleep clust # the same with awake
        end

         colorClusters_all=distinguishable_colors(10);
         A_color=cluster_spatial_footprint_colormap({ntemp},d1,d2,colorClusters_all,group,0.7);
         
         %footprint and region 
         group_rec={};
         CM={};
         for pk=2:10
             [group_rec{1}{pk-1},CM{pk-1},Z,cophList,cophList_adjusted]=cluster_determine_by_suoqin_NMF_firstPeakCoph_stability(ntemp,100,10,pk);
         end
         group_rec{1}{length(unique(group))-1}=group;

%          [~,A_color_region,avg_region,avg_reg_shuf,avg_reg_shuf_95]=DBSCAN_region_quantify_func_simplify_no_plot(group_rec,[2:10],100,length(unique(group))-1,1,{ntemp}); % make sleep clust # the same with awake
         tic;
         A_color={};
         A_color_region={};
         avg_region=[];
         avg_reg_shuf=[];
         avg_reg_shuf_95=[];
         for pk=2:10
%              [A_color{pk},A_color_region{pk},avg_region(pk),avg_reg_shuf(pk),avg_reg_shuf_95(pk)]=DBSCAN_region_quantify_func_single_group_no_plot(group_rec{1}{pk-1},{ntemp},[],{'special'});
            [A_color{pk},A_color_region{pk},avg_region(pk),~,~,~,avg_reg_shuf_temp]=DBSCAN_region_quantify_noShuf_041522(group_rec{1}{pk-1},{ntemp},[]);
            avg_reg_shuf(pk)=mean(avg_reg_shuf_temp);
         end
         toc;
         area_size_summary{i,j}={avg_region,avg_reg_shuf};
         
         clust_idx=length(unique(group))-1;      
%          close
         
         %frame and traj
         subplot(size(behav_vid,2),6,1+(j-1)*6);
         imagesc(frame);
         hold on;
         load(behavExp{i,j});
         plot(behavpos(:,1)*behav.ROI3/behav.trackLength+behav.ROI(1),behavpos(:,2)*behav.ROI3/behav.trackLength+behav.ROI(2),'.','color','c','MarkerSize',5);
         
         % velocity summary
         dist_chg=[];
         for k=1:size(behavpos,1)-1
             dist_chg(k,1)=norm(behavpos(k+1,:)-behavpos(k,:));
         end
         velo_summary_t=dist_chg./(diff(btime(1:end))/1000);%mm/sec
         velo_summary_t(velo_summary_t==inf)=[];
         velo_summary_t(isnan(velo_summary_t))=[];
         velo_summary(i,j)=nanmean(velo_summary_t);
         if j==2
             velo_summary(i,j+1)=nanmedian(velo_summary_t);
         end
         
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

%         individual neuron trace
%          colorClusters2=colorClusters_all;
         group2=group;
%          subplot(size(behav_vid,2),6,3+(j-1)*6);
%          for pk = 1:length(unique(group2))
%             subplot(length(unique(group2)),1,pk);
%             rangee=find(group2 == pk);
%             ctt=1;
%             for k=1:length(rangee)
%                 plot(ntemp.C(rangee(k),1:end)'+(ctt-1)*10,'-','color',colorClusters2(j,:));hold on;
%                 ctt=ctt+1;
%             end
%             axis off;
%          end
         
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
%          plot(avg_reg_shuf_95(2:end),'-^','color',[0.9290    0.6940    0.1250]);
         fr_temp=sum(ntemp.S>0,2)./size(ntemp.S,2);
         fr_temp=fr_temp(fr_temp>10/size(ntemp.S,2)); % remove some small firing rate neurons
         fr_summary{i,j}=mean(fr_temp);
         all_avg_region{i,j}=avg_region;
         
%          [~,~,~,~,~,~,~,intra_all{i,j},inter_all{i,j},intra_shuffle_all{i,j}]=intra_inter_cluster_corr_dis({ntemp},group2,1,'corr');
%          [~,~,~,~,~,~,~,intra_all_dis{i,j},inter_all_dis{i,j},intra_shuffle_all_dis{i,j}]=intra_inter_cluster_corr_dis({ntemp},group2,1,'dis');

    end
end

fr_summary_mat=cell2mat(fr_summary);
%% pairwise correlation, intra/inter
for i=1:3
    corr_aw_slp(i-1,1)=mean(mean(pcorr_summary{i,1}(pcorr_summary{i,1}>0)));
    corr_aw_slp(i-1,2)=mean(mean(pcorr_summary{i,2}(pcorr_summary{i,2}>0)));
end

%% size
for i=1:3
    size_aw_slp(i-1,1)=mean(area_size_summary{i,1}(area_size_summary{i,1}>0));
    size_aw_slp(i-1,2)=mean(area_size_summary{i,2}(area_size_summary{i,2}>0));
end

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
[h,p1]=kstest2(k1(1:10:end),k2(1:50:end));
[h,p2]=kstest2(k1(1:10:end),k3(1:10000:end));
[h,p3]=kstest2(k2(1:50:end),k3(1:10000:end));

mean(cell2mat(intra_all(:,2)))
std(cell2mat(intra_all(:,2)))/length(cell2mat(intra_all(:,2)))^0.5
mean(cell2mat(inter_all(:,2)))
std(cell2mat(inter_all(:,2)))/length(cell2mat(inter_all(:,2)))^0.5
mean(cell2mat(intra_shuffle_all(:,2)))
std(cell2mat(intra_shuffle_all(:,2)))/length(cell2mat(intra_shuffle_all(:,2)))^0.5

k1=cell2mat(intra_all(:,2));
k2=cell2mat(inter_all(:,2));
k3=cell2mat(intra_shuffle_all(:,2));
[h,p1]=kstest2(k1(1:10:end),k2(1:50:end));
[h,p2]=kstest2(k1(1:10:end),k3(1:10000:end));
[h,p3]=kstest2(k2(1:50:end),k3(1:10000:end));

%% intra-inter-cluster distance
intra_all_dis={};
inter_all_dis={};
intra_shuffle_all_dis={};
for i=1:3
    for j=1:size(behav_vid,2)
        load(neuronDat{i,j});
        ntemp=neuron.copy;
        if j==1
            ntemp.C=ntemp.C(:,1:ntemp.num2read(2));
        else
            ntemp.C=ntemp.C(:,ntemp.num2read(2)+1:end);
        end
        ntemp.S=C_to_peakS(ntemp.C);
       
        thresh=3*std(ntemp.S,[],2);
        ntemp.S=thresh_nC_nS(ntemp.S,thresh);

        [intra_all_dis{i,j},inter_all_dis{i,j},~,~,intra_shuffle_all_dis{i,j},~]=intra_inter_cluster_corr_dis({ntemp},group_summary{i,j},1,'dis');

    end
end

subplot(121)

colorClusters_all=distinguishable_colors(10);
h1=cdfplot(cell2mat(intra_all_dis(:,1)));
hold on;
h2=cdfplot(cell2mat(inter_all_dis(:,1)));
h3=cdfplot(cell2mat(intra_shuffle_all_dis(:,1)));
set(h1,'color',[0.4940,0.1840,0.5560]);
set(h2,'color',[0.4660,0.6740,0.1880]);
set(h3,'color',[0.3010,0.7450,0.9330]);

p1=infer_cdf_loc(cell2mat(intra_all_dis(:,1)),mean(cell2mat(intra_all_dis(:,1))));
p2=infer_cdf_loc(cell2mat(inter_all_dis(:,1)),mean(cell2mat(inter_all_dis(:,1))));
p3=infer_cdf_loc(cell2mat(intra_shuffle_all_dis(:,1)),mean(cell2mat(intra_shuffle_all_dis(:,1))));

plot(p1(1),p1(2),'.','color',colorClusters_all(4,:),'MarkerSize',20)
plot(p2(1),p2(2),'.','color',colorClusters_all(5,:),'MarkerSize',20)
plot(p3(1),p3(2),'.','color',colorClusters_all(6,:),'MarkerSize',20)

subplot(122)

h1=cdfplot(cell2mat(intra_all_dis(:,2)));
hold on;
h2=cdfplot(cell2mat(inter_all_dis(:,2)));
h3=cdfplot(cell2mat(intra_shuffle_all_dis(:,2)));
set(h1,'color',[0.4940,0.1840,0.5560]);
set(h2,'color',[0.4660,0.6740,0.1880]);
set(h3,'color',[0.3010,0.7450,0.9330]);

p1=infer_cdf_loc(cell2mat(intra_all_dis(:,2)),mean(cell2mat(intra_all_dis(:,2))));
p2=infer_cdf_loc(cell2mat(inter_all_dis(:,2)),mean(cell2mat(inter_all_dis(:,2))));
p3=infer_cdf_loc(cell2mat(intra_shuffle_all_dis(:,2)),mean(cell2mat(intra_shuffle_all_dis(:,2))));

plot(p1(1),p1(2),'.','color',colorClusters_all(4,:),'MarkerSize',20)
plot(p2(1),p2(2),'.','color',colorClusters_all(5,:),'MarkerSize',20)
plot(p3(1),p3(2),'.','color',colorClusters_all(6,:),'MarkerSize',20)

mean(cell2mat(intra_all_dis(:,1)))
std(cell2mat(intra_all_dis(:,1)))/length(cell2mat(intra_all_dis(:,1)))^0.5
mean(cell2mat(inter_all_dis(:,1)))
std(cell2mat(inter_all_dis(:,1)))/length(cell2mat(inter_all_dis(:,1)))^0.5
mean(cell2mat(intra_shuffle_all_dis(:,1)))
std(cell2mat(intra_shuffle_all_dis(:,1)))/length(cell2mat(intra_shuffle_all_dis(:,1)))^0.5

k1=cell2mat(intra_all_dis(:,1));
k2=cell2mat(inter_all_dis(:,1));
k3=cell2mat(intra_shuffle_all_dis(:,1));
[h,p1]=kstest2(k1(1:10:end),k2(1:50:end));
[h,p2]=kstest2(k1(1:10:end),k3(1:10000:end));
[h,p3]=kstest2(k2(1:50:end),k3(1:10000:end));

mean(cell2mat(intra_all_dis(:,2)))
std(cell2mat(intra_all_dis(:,2)))/length(cell2mat(intra_all_dis(:,2)))^0.5
mean(cell2mat(inter_all_dis(:,2)))
std(cell2mat(inter_all_dis(:,2)))/length(cell2mat(inter_all_dis(:,2)))^0.5
mean(cell2mat(intra_shuffle_all_dis(:,2)))
std(cell2mat(intra_shuffle_all_dis(:,2)))/length(cell2mat(intra_shuffle_all_dis(:,2)))^0.5

k1=cell2mat(intra_all_dis(:,2));
k2=cell2mat(inter_all_dis(:,2));
k3=cell2mat(intra_shuffle_all_dis(:,2));
[h,p1]=kstest2(k1(1:10:end),k2(1:50:end));
[h,p2]=kstest2(k1(1:10:end),k3(1:10000:end));
[h,p3]=kstest2(k2(1:50:end),k3(1:10000:end));

%% region size compare, awake/sleep
awake_ori_size=[];
awake_shuf_size=[];
sleep_ori_size=[];
sleep_shuf_size=[];

for i=1:3
    awake_ori_size=[awake_ori_size;original_area_size_summary{i,1}{1}];
    awake_shuf_size=[awake_shuf_size;nanmean(original_area_size_summary{i,1}{2})];
    sleep_ori_size=[sleep_ori_size;original_area_size_summary{i,2}{1}];
    sleep_shuf_size=[sleep_shuf_size;nanmean(original_area_size_summary{i,2}{2})];    
end