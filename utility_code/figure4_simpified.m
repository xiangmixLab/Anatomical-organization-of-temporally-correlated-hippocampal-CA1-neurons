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

        [group,CMO,Z,cophList,cophList_adjusted]=cluster_determine_by_suoqin_NMF_firstPeakCoph_stability(ntemp,100,10,[]);
        [A_color,A_color_region,original_area_size_summary{i,j}{1},~,~,~,original_area_size_summary{i,j}{2},~]=DBSCAN_region_quantify_022422(group,{ntemp},[]); % make sleep clust # the same with awake      

         %colorClusters_all=distinguishable_colors(10);
         %A_color=cluster_spatial_footprint_colormap({ntemp},d1,d2,colorClusters_all,group,0.7);

         clust_idx=length(unique(group))-1;      

         %frame and traj
         subplot(size(behav_vid,2),4,1+(j-1)*4);
         imagesc(frame);
         hold on;
         load(behavExp{i,j});
         plot(behavpos(:,1)*behav.ROI3/behav.trackLength+behav.ROI(1),behavpos(:,2)*behav.ROI3/behav.trackLength+behav.ROI(2),'color','c','MarkerSize',5);
        

         %individual neuron trace
%          subplot(size(behav_vid,2),4,2+(j-1)*4);
%          group2=group;
%          clustered_timeSeries2_panel(ntemp.C,group2)

         % raster map
         subplot(size(behav_vid,2),4,2+(j-1)*4);
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
         caxis([0 80]);
         hold on
         for p = 2:length(positionC2)-1
             line(get(gca,'XLim'),[positionC2(p) positionC2(p)],'LineWidth',1,'Color','r'); hold on;
         end
         set(gca,'FontSize',8)
         labels = strcat('C',cellstr(num2str([1:optimalK2]')));
         ytick0 = positionC2(1:end-1)+diff(positionC2)/2;
         set(gca,'Ytick',ytick0);set(gca,'YtickLabel',labels,'FontName','Arial','FontSize',10)
         xlabel('Frames','FontSize',10)
         ylabel('Neurons','FontSize',10)

         %footprint
         subplot(size(behav_vid,2),4,3+(j-1)*4);
         imagesc(A_color);
         
         %region
         subplot(size(behav_vid,2),4,4+(j-1)*4);
         imagesc(A_color_region);
        
         
%          [~,~,~,~,~,~,~,intra_all{i,j},inter_all{i,j},intra_shuffle_all{i,j}]=intra_inter_cluster_corr_dis({ntemp},group2,1,'corr');
%          [~,~,~,~,~,~,~,intra_all_dis{i,j},inter_all_dis{i,j},intra_shuffle_all_dis{i,j}]=intra_inter_cluster_corr_dis({ntemp},group2,1,'dis');

    end
end
