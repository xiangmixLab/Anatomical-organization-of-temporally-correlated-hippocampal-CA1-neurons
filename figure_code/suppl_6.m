%% Supplemental Fig6: supporting data for Fig2 and suppl5

load('D:\final_HDAC_AD_automatic_processing\dataloc\Yanjun_nnrevision_circle_square_102919.mat')
%% check group numbers
group_all=[];
for i=1:length(foldername)
    load([foldername{i},'\','square_placement_cell_cluster','\','variables_clustering_0.35PC.mat']);
    group_all(i,1)=length(unique(group{end}));
    load([foldername{i},'\','circle_placement_cell_cluster','\','variables_clustering_0.35PC.mat']);
    group_all(i,2)=length(unique(group{end}));
end

%% panel a: behav trajectory
load([foldername{1},'\','square_Behav.mat']);
plot(behav.position(:,1),behav.position(:,2),'color',[0.5,0.5,0.5],'lineWidth',2);


%% Panel b: cluster ensemble, region, size change wit h# of clusters
% M3422F M3423F M3426F
mice_all=[1 8 11];
load('C:\Users\exx\Desktop\HDAC paper fig and method\SFN2019\SFN2019 fig\SCI clustNum plot\cir box\fig1_dataset.mat');
group_record(2,:)=group_record;

for i=1:3
    mice=mice_all(i);
    clust_num=[2:10];
    shuffles_num=100;
    gp_rec=group_record(:,mice);
    
    load([foldername{mice},'\','neuronIndividuals_new.mat'])
    load([foldername{mice},'\','square_placement_cell_cluster','\','variables_clustering_0.35PC.mat'])
    group_model=group{1};
    
    clust_idx=length(unique(group_model));
    cond_sign=1;

    gp_rec{cond_sign}{clust_idx-1}=group_model;

    DBSCAN_region_quantify_func_simplify(gp_rec,clust_num,shuffles_num,clust_idx-1,cond_sign,neuronIndividuals_new)
end
%% panel c: intra, intra shuffle and inter clust corr, all mice
for tk=1:12
    cond_sign=1;
    load([foldername{tk},'\','neuronIndividuals_new.mat'])
    load([foldername{tk},'\','square_placement_cell_cluster','\','variables_clustering_0.35PC.mat'])
    group_model=group{cond_sign};
    [intra_all{tk},inter_all{tk},~,~,intra_shuffle_all{tk},~]=intra_inter_cluster_corr_dis(neuronIndividuals_new,group_model,cond_sign,'corr');
end

colorClusters_all=colormap(lines);close

a1x=figure;
% subplot(211)

h1=cdfplot(cell2mat(intra_all'));
hold on;
h2=cdfplot(cell2mat(inter_all'));
h3=cdfplot(cell2mat(intra_shuffle_all'));
set(h1,'color',[0.4940,0.1840,0.5560]);
set(h2,'color',[0.4660,0.6740,0.1880]);
set(h3,'color',[0.3010,0.7450,0.9330]);

p1=infer_cdf_loc(cell2mat(intra_all'),mean(cell2mat(intra_all')));
p2=infer_cdf_loc(cell2mat(inter_all'),mean(cell2mat(inter_all')));
p3=infer_cdf_loc(cell2mat(intra_shuffle_all'),mean(cell2mat(intra_shuffle_all')));

plot(p1(1),p1(2),'.','color',colorClusters_all(4,:),'MarkerSize',20)
plot(p2(1),p2(2),'.','color',colorClusters_all(5,:),'MarkerSize',20)
plot(p3(1),p3(2),'.','color',colorClusters_all(6,:),'MarkerSize',20)

mean(cell2mat(intra_all'))
std(cell2mat(intra_all'))/length(cell2mat(intra_all')).^0.5
mean(cell2mat(inter_all'))
std(cell2mat(inter_all'))/length(cell2mat(inter_all')).^0.5
mean(cell2mat(intra_shuffle_all'))
std(cell2mat(intra_shuffle_all'))/length(cell2mat(intra_shuffle_all')).^0.5


k1=cell2mat(intra_all');
k2=cell2mat(inter_all');
k3=cell2mat(intra_shuffle_all');

[h,p1]=kstest2(k1(1:100:end),k2(1:100:end)) % confirmed that k1 and k1(1:100:end) can be considered coming from the same distribution using kstest2, and their shape is very similar
[h,p2]=kstest2(k1(1:100:end),k3(1:10000:end))
%% panel d: distance-temporal correlation plot, all mice
mice_all=[1:12];
cluster_file_name={};
for tk=1:12
    fname=foldername{mice_all(tk)};
    cluster_file_name{tk}=[fname,'\','square_placement_cell_cluster','\','variables_clustering_0.35PC.mat'];
    session=1;
end

[all_dis_corr_circle_all_intra,all_dis_corr_circle_all_inter,all_dis_corr_circle_all_all,midx,f,gof]=pairwise_dis_tempCorr_091421(foldername,cluster_file_name,session);

[H, pValue, KSstatistic] = kstest_2s_2d(all_dis_corr_circle_all_intra(1:1000:end,:), all_dis_corr_circle_all_inter(1:2000:end,:));

mean(all_dis_corr_circle_all_intra(:,1))
std(all_dis_corr_circle_all_intra(:,1))/length(all_dis_corr_circle_all_intra(:,1))^0.5

mean(all_dis_corr_circle_all_inter(:,1))
std(all_dis_corr_circle_all_inter(:,1))/length(all_dis_corr_circle_all_inter(:,1))^0.5

mean(all_dis_corr_circle_all_all(:,1))
std(all_dis_corr_circle_all_all(:,1))/length(all_dis_corr_circle_all_all(:,1))^0.5


[h,p]=kstest2(all_dis_corr_circle_all_all(1:100:end,1),all_dis_corr_circle_all_all(1:100:end,1))

[rho,pval] = corr(all_dis_corr_circle_all_all(1:100:end,1),all_dis_corr_circle_all_all(1:100:end,2),'Type','spearman','Rows','complete')

%% panel e: cluster region size change with cluster number, all mice
avg_max_region={};
avg_max_reg_shuf_all={};
avg_max_reg_shuf_all_95={};

% group_record={};
% for tk=1:12
%     for j=1:2
%         for k=2:10
%             load([foldername{tk},'\','neuronIndividuals_new.mat']);
%             [group,CM,Z,cophList,cophList_adjusted]=cluster_determine_by_suoqin_NMF_firstPeakCoph_stability(neuronIndividuals_new{j},100,10,k);
%             group_record{j,tk}{k}=group;
%         end 
%     end
% end

load('C:\Users\exx\Desktop\HDAC paper fig and method\SFN2019\SFN2019 fig\SCI clustNum plot\cir box\fig1_dataset_group.mat');

for tk=1:12
    load([foldername{tk},'\neuronIndividuals_new.mat']);
    gp_rec=group_record(:,tk)
    clust_num=[2:10];
    shuffles_num=100;
    cond_sign=1;

    load([foldername{tk},'\square_placement_cell_cluster','\','variables_clustering_0.35PC.mat'])
    group_model=group{1};
    clust_idx=length(unique(group_model));
    gp_rec{cond_sign}{clust_idx-1}=group_model;

    [~,~,avg_max_region{tk},avg_max_reg_shuf_all{tk},avg_max_reg_shuf_all_95{tk}]=DBSCAN_region_quantify_func_simplify(gp_rec,clust_num,shuffles_num,clust_idx-1,cond_sign,neuronIndividuals_new);
    close
end
save('C:\Users\exx\Desktop\HDAC paper fig and method\SFN2019\SFN2019 fig\Round10\fig2 square panels\d_data.mat','avg_max_region','avg_max_reg_shuf_all','avg_max_reg_shuf_all_95')

for tk=1:12
   avg_max_region_mat(tk,:)=avg_max_region{tk}(1,:);
   avg_max_reg_shuf_all_mat(tk,:)=avg_max_reg_shuf_all{tk}(1,:);
   avg_max_reg_shuf_all_95_mat(tk,:)=avg_max_reg_shuf_all_95{tk}(1,:);
end

mean_avg_max_region_mat=mean(avg_max_region_mat,1);
se_avg_max_region_mat=std(avg_max_region_mat,[],1)/(size(avg_max_region_mat,1)^0.5);
min_avg_max_region_mat=mean(avg_max_region_mat,1)-se_avg_max_region_mat;
max_avg_max_region_mat=mean(avg_max_region_mat,1)+se_avg_max_region_mat;

mean_avg_max_reg_shuf_all_mat=mean(avg_max_reg_shuf_all_mat,1);
se_avg_max_reg_shuf_all_mat=std(avg_max_reg_shuf_all_mat,[],1)/(size(avg_max_reg_shuf_all_mat,1)^0.5);
min_avg_max_reg_shuf_all_mat=mean(avg_max_reg_shuf_all_mat,1)-se_avg_max_reg_shuf_all_mat;
max_avg_max_reg_shuf_all_mat=mean(avg_max_reg_shuf_all_mat,1)+se_avg_max_reg_shuf_all_mat;

mean_avg_max_reg_shuf_all_95_mat=mean(avg_max_reg_shuf_all_95_mat,1);
se_avg_max_reg_shuf_all_95_mat=std(avg_max_reg_shuf_all_95_mat,[],1)/(size(avg_max_reg_shuf_all_95_mat,1)^0.5);
min_avg_max_reg_shuf_all_95_mat=mean(avg_max_reg_shuf_all_95_mat,1)-se_avg_max_reg_shuf_all_95_mat;
max_avg_max_reg_shuf_all_95_mat=mean(avg_max_reg_shuf_all_95_mat,1)+se_avg_max_reg_shuf_all_95_mat;

x=[[0:length(mean_avg_max_region_mat)-1]',[0:length(mean_avg_max_reg_shuf_all_mat)-1]',[0:length(mean_avg_max_reg_shuf_all_95_mat)-1]'];
y=[mean_avg_max_region_mat',mean_avg_max_reg_shuf_all_mat',mean_avg_max_reg_shuf_all_95_mat'];
minMaxRange={[min_avg_max_region_mat',max_avg_max_region_mat'],[min_avg_max_reg_shuf_all_mat',max_avg_max_reg_shuf_all_mat'],[min_avg_max_reg_shuf_all_95_mat',max_avg_max_reg_shuf_all_95_mat']};
% plotwithMaxMixRange(x,y,minMaxRange,{[0 0 0.9],[0 0 0.55],[0 0 0.20]},{'-o','-s','-^'});

errorbar(0:length(mean_avg_max_region_mat)-1,mean_avg_max_region_mat,se_avg_max_region_mat,'-o','color',[0    0.4470    0.7410]);
hold on;
errorbar(0:length(mean_avg_max_reg_shuf_all_mat)-1,mean_avg_max_reg_shuf_all_mat,se_avg_max_reg_shuf_all_mat,'-s','color',[0.8500    0.3250    0.0980]);
errorbar(0:length(mean_avg_max_reg_shuf_all_95_mat)-1,mean_avg_max_reg_shuf_all_95_mat,se_avg_max_reg_shuf_all_95_mat,'-^','color',[0.9290    0.6940    0.1250]);

xlim([-1 9])

%% panel f: ensemble map
for tk=1:3
    mice=mice_all(tk);
    load([foldername{mice},'\','square_placement_cell_cluster','\','variables_clustering_0.35PC.mat'])
    group=group{1};       
    uni_group=unique(group);
    pct=1;
    figure;
    for i=1:length(unique(group))

        subplot(4,length(unique(group)),pct);
        hold on;
        colormap(jet)

        load([foldername{mice},'\','square_placement_cell_cluster','\','cluster_ensemble_analysis','\','cluster',num2str(i),'_neuron_comparingfiringRate_S_binsize10data_S.mat'])
        t1=double(firingRateSmoothing2);
        t1(countTime==0)=nan;
        pcolor(t1);
        colorbar

        idxx=find(group==uni_group(i));
        load([foldername{mice},'\','square_results\single_cell_firing_profile_S.mat'])
        max_inds={};
        peak_rates={};
        PF_radius={};
        for j=1:1:length(idxx)
            t1=firingrateS{idxx(j)};
            t1(countTime==0)=nan;
            [max_inds{j},peak_rates{j},PF_radius{j}]=find_max_response_loc(t1);  
            max_peak_rates(j)=max(peak_rates{j});
        end

        max_max_peak_rates=max(max_peak_rates);
        for j=1:1:length(idxx)
            for k1=1:size(max_inds{j},1)
                cen=max_inds{j}(k1,:);
                scatter1=scatter(round(cen(1,2)),round(cen(1,1)),50,[0.1 0.1 0.1],'.');
            end
        end

        axis ij
        shading flat;
        axis image
        axis off
        
        [cell_fields]=single_field_cell_clust_determine(group,i,firingrateS,countTime);
        
        subplot(4,length(unique(group)),pct+length(unique(group)));
        pcolor(cell_fields{1});
        colorbar
        axis ij
        shading flat;
        axis image
        axis off
        colormap(jet)
        
        subplot(4,length(unique(group)),pct+length(unique(group))*2);
        pcolor(cell_fields{2});
        colorbar
        axis ij
        shading flat;
        axis image
        axis off
        colormap(jet)
        
        subplot(4,length(unique(group)),pct+length(unique(group))*3);
        pcolor(cell_fields{3});
        colorbar
        axis ij
        shading flat;
        axis image
        axis off
        colormap(jet)
        
        pct=pct+1;
    end
    set(gcf,'renderer','painters');
end



% %% Panel g distance-temporal correlation plot
% for tk=1:3
%     fname=foldername{mice};
%     cluster_file_name=[foldername{mice},'\','square_placement_cell_cluster','\','variables_clustering_0.35PC.mat'];
%     [all_dis_corr_circle_all_intra,all_dis_corr_circle_all_inter]=pairwise_dis_tempCorr({fname},cluster_file_name,1);
% end
%% Panel f infoscore, distribution
infoPerSpike={};
for tk=3:3
    figure;
    mice=mice_all(tk);
    fname=foldername{mice};
    cluster_file_name=[foldername{mice},'\','square_placement_cell_cluster','\','variables_clustering_0.35PC.mat'];
   
    neuronName=[foldername{mice},'\','neuronIndividuals_new.mat'];
    behavName=[foldername{mice},'\','square_results\current_condition_behav.mat'];
    threshName=[foldername{mice},'\','thresh_and_ROI.mat'];
    load(cluster_file_name);
    load(neuronName);
    load(behavName);
    load(threshName);
     group1=group{1}
    [firingrateAll,countAll,~,countTime] = calculatingCellSpatialForSingleData_Suoqin(neuronIndividuals_new{1},behavpos,behavtime,maxbehavROI,10,1:size(neuronIndividuals_new{1}.S,1),3*std(neuronIndividuals_new{1}.S,[],2),'S',[],[],[0.1 1000000],10);
    for j=1:length(unique(group1))
        [infoPerSecond{j}, infoPerSpike{j}] = comparisonSpatialInfo_adapt(firingrateAll(group1==j), countAll(group1==j), countTime,[0.1],[],10);
        cdfplot(infoPerSpike{j});
        hold on
    end

    legend({'cluster1','cluster2','cluster3','cluster4','cluster5','cluster6','cluster7','cluster8'});


end

%% new panel g: percentage of place cell
figure;
for tk=1:3  
    mice=mice_all(tk);
    fname=foldername{mice};
    cluster_file_name=[foldername{mice},'\','square_placement_cell_cluster','\','variables_clustering_0.35PC.mat'];
    pc_file_name=[foldername{mice},'\','square_info_score_placement_cell','\','place_cells_info_square_info_scor_binsize10_S_spk.mat'];
  
    neuronName=[foldername{mice},'\','neuronIndividuals_new.mat'];
    behavName=[foldername{mice},'\','square_results\current_condition_behav.mat'];
    threshName=[foldername{mice},'\','thresh_and_ROI.mat'];
    load(neuronName);
    load(cluster_file_name);
    load(pc_file_name);
    group1=group{1};
    
    colorClusters_all=distinguishable_colors(10);
    colorClusters_all(11,:)=[0.8 0.8 0.8];
    
    group2=ones(length(group1),1)*11;
    group2(place_cells)=group1(place_cells);
    
    A_color=cluster_spatial_footprint_colormap({neuronIndividuals_new{2}},240,376,colorClusters_all,group2);
    
    percentage_pc=[];
    num_pc=[];
    infoo=[];
    pc_idx=zeros(length(group1),1);
    pc_idx(place_cells)=1;
    for i=1:length(unique(group1))
        idx=(group1==i).*pc_idx;
        percentage_pc(i)=sum(idx)/sum(group1==i);
        num_pc(i)=sum(idx);
        infoo(i)=mean(TinfoPerSecond.infoScore(logical(idx)));
    end

    subplot(length(mice_all),4,4*tk-3)
    imagesc(A_color);
    subplot(length(mice_all),4,4*tk-2)
    b1=bar(percentage_pc','FaceColor','flat');
    for i=1:length(percentage_pc)
        b1.CData(i,:)=[colorClusters_all(i,:)]
    end
    subplot(length(mice_all),4,4*tk-1)
    b2=bar(num_pc,'FaceColor','flat');
    for i=1:length(percentage_pc')
        b2.CData(i,:)=[colorClusters_all(i,:)]
    end
    subplot(length(mice_all),4,4*tk)
    b3=bar(infoo,'FaceColor','flat');
    for i=1:length(percentage_pc')
        b3.CData(i,:)=[colorClusters_all(i,:)]
    end
end

%% panel g new: pairwise temporal corr-pairwise rate map correlation (pearson)
figure;
for tk=1:3
    mice=mice_all(tk);
    load([foldername{mice},'\neuronIndividuals_new.mat']);
    load([foldername{mice},'\','square_results','\','single_cell_firing_profile_S.mat']);
    temp_corr=1-pdist(neuronIndividuals_new{1}.C,'correlation');
    rateMap_corr=[];
    ctt=1;
    idx=zeros(length(temp_corr),1);
    for i=1:length(firingrateS)-1
        for j=i+1:length(firingrateS)
            cell1=firingrateS{i};
            cell2=firingrateS{j};
            if ~isempty(cell1)&&~isempty(cell2)
                cell1=filter2DMatrices(cell1,1);
                cell2=filter2DMatrices(cell2,1);
                cell1(countTime==0)=-inf; % get rid of non trespassing periods, avoid unnecessary correlation boost
                cell2(countTime==0)=-inf;
                cell1=reshape(cell1,size(cell1,1)*size(cell1,2),1);
                cell2=reshape(cell2,size(cell2,1)*size(cell2,2),1);
                cell1(cell1==-inf)=[]; % get rid of non trespassing periods, avoid unnecessary correlation boost
                cell2(cell2==-inf)=[];

                rm_corr=corrcoef(cell1,cell2);
                rateMap_corr(ctt)=rm_corr(1,2);                
            else
                idx(ctt)=1;
            end
            ctt=ctt+1;
        end
     end
     

%      f1=fit(rateMap_corr',temp_corr','power2');
     subplot(1,3,tk)
     rateMap_corr=rateMap_corr(~logical(idx));
     temp_corr=temp_corr(~logical(idx));
     h1=plot(rateMap_corr(1:50:end),temp_corr(1:50:end),'.');
     set([h1],'color',[0.5 0.5 0.5],'MarkerSize',10)
end

%% new d: pairwise temporal corr-pairwise rate map correlation (pearson)
figure;
gof_collection={};
for tk=1:3
    mice=mice_all(tk);
    load([foldername{mice},'\neuronIndividuals_new.mat']);
    load([foldername{mice},'\','square_results','\','single_cell_firing_profile_S.mat']);
    temp_corr=1-pdist(neuronIndividuals_new{1}.C,'correlation');
    rateMap_corr=[];
    ctt=1;
    for i=1:length(firingrateS)-1
        for j=i+1:length(firingrateS)
            cell1=firingrateS{i};
            cell2=firingrateS{j};
            cell1=filter2DMatrices(cell1,1);
            cell2=filter2DMatrices(cell2,1);
            fr(countTime==0)=-inf; % get rid of non trespassing periods, avoid unnecessary correlation boost

            cell1=reshape(cell1,size(cell1,1)*size(cell1,2),1);
            cell2=reshape(cell2,size(cell2,1)*size(cell2,2),1);
            if ~isempty(cell1)&&~isempty(cell2)
                rm_corr=corrcoef(cell1,cell2);
                rateMap_corr(ctt)=rm_corr(1,2);
                ctt=ctt+1;
            else
                rateMap_corr(ctt)=nan;
                ctt=ctt+1;
            end
        end
     end
     

%      f1=fit(rateMap_corr',temp_corr','power2');
     subplot(1,3,tk)
     nan_idx=(isnan(rateMap_corr')+isnan(temp_corr'))>0;
     rateMap_corr=rateMap_corr(~nan_idx);
     temp_corr=temp_corr(~nan_idx);
     [f1,gof]=fit(rateMap_corr',temp_corr','poly1');
     h1=plot(f1,rateMap_corr(1,1:40:end)',temp_corr(1,1:40:end)');
     set([h1],'color',[0.5 0.5 0.5],'MarkerSize',10)
     gof_collection{tk,1}=gof;
     [rho,p]=corr(rateMap_corr(1:100:end)',temp_corr(1:100:end)');
     gof_collection{tk,2}=[rho,p]
end

