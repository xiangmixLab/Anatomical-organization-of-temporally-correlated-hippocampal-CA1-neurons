% 
run('D:\Xu_clusterting_paper_prep11_2020\final code\data_prepare\data_info.m')
all_colors=distinguishable_colors(10);
load(['D:\Xu_clusterting_paper_prep11_2020\final data\final_cluster_data\Fig2_cir_rec_clust_original.mat'])

mice=[1,8,11];

%% A: cluster footprint
for i=1:length(mice)
    load([foldername_fig2{mice(i)},'\','neuronIndividuals_new.mat'])
    [A_color1,A_color_region1]=DBSCAN_region_quantify_022422(group_ori_fig2{mice(i),2},neuronIndividuals_new,[]);
    
    subplot(3,2,2*i-1);
    imagesc(A_color1);
    subplot(3,2,2*i);
    imagesc(A_color_region1);
end

%% panel e: ensemble map
mice_all=[1,8,11];
peak_rate={};
for tk=1:3
    mice=mice_all(tk);
    load([foldername{mice},'\','circle_placement_cell_cluster','\','variables_clustering_0.35PC.mat'])
    group=group{2};       
    uni_group=unique(group);
    pct=1;
    figure;
    pkrate_temp=[];
    for i=1:length(unique(group))
        
        % ensemble ratemap
        subplot(5,length(unique(group)),pct);
        hold on;
        colormap(jet)

        load([foldername{mice},'\','circle_placement_cell_cluster','\','cluster_ensemble_analysis','\','cluster',num2str(i),'_neuron_comparingfiringRate_S_binsize10data_S.mat'])
        t1=double(firingRateSmoothing2);
        t1(countTime==0)=nan;
        pcolor(t1);
        colorbar
        axis ij
        shading flat;
        axis image
        axis off
        pkrate_temp(i)=max(firingRateSmoothing2(:));
        % peak ratemap pos
        subplot(5,length(unique(group)),pct+length(unique(group)));hold on;
        idxx=find(group==uni_group(i));
        load([foldername{mice},'\','circle_results\single_cell_firing_profile_S.mat'])

        cen=[];
        for j=1:1:length(idxx)
            tt=firingrateS{idxx(j)};
            tt=filter2DMatrices(tt,1);

            [cen(j,1),cen(j,2)]=find(tt==max(tt(:)));
        end
        scatter(round(cen(:,2)),round(cen(:,1)),50,[0.1 0.1 0.1],'.','MarkerFaceAlpha',.5);
        axis ij
        shading flat;
        axis image
        axis off
       
        
        % single cells
        
        pct=pct+1;
    end
    peak_rate{tk}=pkrate_temp;
    set(gcf,'renderer','painters');
end

%% panel a: behav trajectory
load([foldername_fig2{1},'\','circle_Behav.mat']);
plot(behav.position(:,1),behav.position(:,2),'color',[0.5,0.5,0.5],'lineWidth',2);

%% Panel b: cluster ensemble, region, size change wit h# of clusters
% M3422F M3423F M3426F
mice_all=[1 8 11];
load('C:\Users\exx\Desktop\HDAC paper fig and method\SFN2019\SFN2019 fig\SCI clustNum plot\cir box\fig1_dataset.mat');
group_record(2,:)=group_record;

for i=1:12
    mice=mice_all(i);
    clust_num=[2:10];
    shuffles_num=100;
    gp_rec=group_record(:,mice);
    
    load([foldername_fig2{mice},'\','neuronIndividuals_new.mat'])
    load([foldername_fig2{mice},'\','circle_placement_cell_cluster','\','variables_clustering_0.35PC.mat'])
    group_model=group{2};
    
    clust_idx=length(unique(group_model));
    cond_sign=2;

    gp_rec{cond_sign}{clust_idx-1}=group_model;

    DBSCAN_region_quantify_func_simplify(gp_rec,clust_num,shuffles_num,clust_idx-1,cond_sign,neuronIndividuals_new)
end
%% panel c: intra, inter, intra shuffle pdf, all mice
for tk=1:12
    cond_sign=2;
    load([foldername_fig2{tk},'\','neuronIndividuals_new.mat'])
    load([foldername_fig2{tk},'\','circle_placement_cell_cluster','\','variables_clustering_0.35PC.mat'])
    group_model=group{2};
    [~,~,~,~,~,~,~,intra_all{tk},inter_all{tk},intra_shuffle_all{tk}]=intra_inter_cluster_corr_dis(neuronIndividuals_new,group_model,cond_sign,'corr');
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

k1=cell2mat(intra_all');
k2=cell2mat(inter_all');
k3=cell2mat(intra_shuffle_all');

[h,p1]=kstest2(k1(1:100:end),k2(1:100:end))
[h,p2]=kstest2(k1(1:100:end),k3(1:10000:end))

%% panel d: cluster region size change with cluster number, all mice
avg_max_region={};
avg_max_reg_shuf_all={};
avg_max_reg_shuf_all_95={};

group_record={};
for tk=1:12
    for j=1:2
        for k=2:10
            load([foldername_fig2{tk},'\','neuronIndividuals_new.mat']);
            [group,CM,Z,cophList,cophList_adjusted]=cluster_determine_by_suoqin_NMF_firstPeakCoph_stability(neuronIndividuals_new{j},100,10,k);
            group_record{j,tk}{k}=group;
        end 
    end
end

load('C:\Users\exx\Desktop\HDAC paper fig and method\SFN2019\SFN2019 fig\SCI clustNum plot\cir box\fig1_dataset_group.mat');

for tk=1:12
    load([foldername_fig2{tk},'\neuronIndividuals_new.mat']);
    gp_rec=group_record(:,tk)
    clust_num=[2:10];
    shuffles_num=100;
    cond_sign=2;

    load([foldername_fig2{tk},'\circle_placement_cell_cluster','\','variables_clustering_0.35PC.mat'])
    group_model=group{2};
    clust_idx=length(unique(group_model));
    gp_rec{cond_sign}{clust_idx-1}=group_model;

    [~,~,avg_max_region{tk},avg_max_reg_shuf_all{tk},avg_max_reg_shuf_all_95{tk}]=DBSCAN_region_quantify_func_simplify(gp_rec,clust_num,shuffles_num,clust_idx-1,cond_sign,neuronIndividuals_new);
    close
end

save('C:\Users\exx\Desktop\HDAC paper fig and method\SFN2019\SFN2019 fig\Round10\fig2 panels\d_data.mat','avg_max_region','avg_max_reg_shuf_all','avg_max_reg_shuf_all_95')
for tk=1:12
   avg_max_region_mat(tk,:)=avg_max_region{tk}(2,:);
   avg_max_reg_shuf_all_mat(tk,:)=avg_max_reg_shuf_all{tk}(2,:);
   avg_max_reg_shuf_all_95_mat(tk,:)=avg_max_reg_shuf_all_95{tk}(2,:);
end
%load('C:\Users\exx\Desktop\HDAC paper fig and method\SFN2019\SFN2019 fig\Round9\fig2 panels\d.mat')
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

%% supplemental panel d: pdis, pcorr relationship
% mice_all=[1:12];
% cluster_file_name={};
% for tk=1:12
%     fname=foldername_fig2{mice_all(tk)};
%     load([foldername_fig2{tk},'\','neuronIndividuals_new.mat'])
%     ntemp_all{tk}=neuronIndividuals_new{2}.copy;
%     cluster_file_name{tk}=[fname,'\','circle_placement_cell_cluster','\','variables_clustering_0.35PC.mat'];
%     load(cluster_file_name{tk});
%     group_all{tk}=group{2};
% end

session=2;
[all_dis_corr_circle_all_intra,all_dis_corr_circle_all_inter,all_dis_corr_circle_all_all,midx,f,gof]=pairwise_dis_tempCorr_091421(foldername_fig2,clust_filename(:,session),session);

[H, pValue, KSstatistic] = kstest_2s_2d(all_dis_corr_circle_all_intra(1:1000:end,:), all_dis_corr_circle_all_inter(1:2000:end,:));

mean(all_dis_corr_circle_all_intra(:,1))
std(all_dis_corr_circle_all_intra(:,1))/length(all_dis_corr_circle_all_intra(:,1))^0.5

mean(all_dis_corr_circle_all_inter(:,1))
std(all_dis_corr_circle_all_inter(:,1))/length(all_dis_corr_circle_all_inter(:,1))^0.5

mean(all_dis_corr_circle_all_all(:,1))
std(all_dis_corr_circle_all_all(:,1))/length(all_dis_corr_circle_all_all(:,1))^0.5


[h,p]=kstest2(all_dis_corr_circle_all_all(1:100:end,1),all_dis_corr_circle_all_all(1:100:end,1))

[rho,pval] = corr(all_dis_corr_circle_all_all(1:100:end,1),all_dis_corr_circle_all_all(1:100:end,2),'Type','spearman','Rows','complete')
%% panel e: ensemble map
peak_rate={};
for tk=1:3
    mice=mice_all(tk);
    load([foldername_fig2{mice},'\','circle_placement_cell_cluster','\','variables_clustering_0.35PC.mat'])
    group=group{2};       
    uni_group=unique(group);
    pct=1;
    figure;
    for i=1:length(unique(group))

        subplot(4,length(unique(group)),pct);
        hold on;
        colormap(jet)

        load([foldername_fig2{mice},'\','circle_placement_cell_cluster','\','cluster_ensemble_analysis','\','cluster',num2str(i),'_neuron_comparingfiringRate_S_binsize10data_S.mat'])
        t1=double(firingRateSmoothing2);
        t1(countTime==0)=nan;
        pcolor(t1);
        colorbar

        idxx=find(group==uni_group(i));
        load([foldername_fig2{mice},'\','circle_results\single_cell_firing_profile_S.mat'])
        max_inds={};
        prates={};
        PF_radius={};
        for j=1:1:length(idxx)
            tt=firingrateS{idxx(j)};
            tt(countTime==0)=nan;
            [max_inds{j},prates{j},PF_radius{j}]=find_max_response_loc(tt);  % find all local peak locations, good for multi-field cells
        end

%         max_max_peak_rates=max(max_peak_rates);
        for j=1:1:length(idxx)
            for k1=1:size(max_inds{j},1)
                cen=max_inds{j}(k1,:);
                scatter1=scatter(round(cen(1,2)),round(cen(1,1)),50,[0.1 0.1 0.1],'.','MarkerFaceAlpha',.5);
            end
        end

        axis ij
        shading flat;
        axis image
        axis off
        
        peak_rate{tk}(1,i)=nanmax(t1(:));
        
        [cell_fields]=single_field_cell_clust_determine(group,i,firingrateS,countTime);
        
        subplot(4,length(unique(group)),pct+length(unique(group)));
        pcolor(cell_fields{1});
        colorbar
        axis ij
        shading flat;
        axis image
        axis off
        colormap(jet)
        peak_rate{tk}(2,i)=nanmax(cell_fields{1}(:));
         
        subplot(4,length(unique(group)),pct+length(unique(group))*2);
        pcolor(cell_fields{2});
        colorbar
        axis ij
        shading flat;
        axis image
        axis off
        colormap(jet)
        peak_rate{tk}(3,i)=nanmax(cell_fields{2}(:));
        
        subplot(4,length(unique(group)),pct+length(unique(group))*3);
        pcolor(cell_fields{3});
        colorbar
        axis ij
        shading flat;
        axis image
        axis off
        colormap(jet)
        peak_rate{tk}(4,i)=nanmax(cell_fields{3}(:));
        
        pct=pct+1;
    end
    set(gcf,'renderer','painters');
end



% %% Panel e distance-temporal correlation plot
% for tk=1:3
%     fname=foldername_fig2{mice};
%     cluster_file_name=[foldername_fig2{mice},'\','circle_placement_cell_cluster','\','variables_clustering_0.35PC.mat'];
%     [all_dis_corr_circle_all_intra,all_dis_corr_circle_all_inter]=pairwise_dis_tempCorr({fname},cluster_file_name,2);
% end
%% Panel f infoscore, distribution
for tk=1:1
    figure;
    mice=mice_all(tk);
    fname=foldername_fig2{mice};
    cluster_file_name=[foldername_fig2{mice},'\','circle_placement_cell_cluster','\','variables_clustering_0.35PC.mat'];
   
    neuronName=[foldername_fig2{mice},'\','neuronIndividuals_new.mat'];
    behavName=[foldername_fig2{mice},'\','circle_results\current_condition_behav.mat'];
    threshName=[foldername_fig2{mice},'\','thresh_and_ROI.mat'];
    load(cluster_file_name);
    load(neuronName);
    load(behavName);
    load(threshName);
     group1=group{2}
    [firingrateAll,countAll,~,countTime] = calculatingCellSpatialForSingleData_Suoqin(neuronIndividuals_new{2},behavpos,behavtime,maxbehavROI,10,1:size(neuronIndividuals_new{2}.S,1),3*std(neuronIndividuals_new{2}.S,[],2),'S',[],[],[0.1 1000000],10);
    for j=1:length(unique(group1))
        [infoPerSecond{j}, infoPerSpike{j}] = comparisonSpatialInfo_adapt(firingrateAll(group1==j), countAll(group1==j), countTime,[0.1],[],10);
        cdfplot(infoPerSpike{j});
        hold on
    end

    legend({'cluster1','cluster2','cluster3','cluster4','cluster5'});


end

%% panel g: pairwise temporal corr-pairwise rate map correlation (pearson)
figure;
gof_collection={};
for tk=1:3
    mice=mice_all(tk);
    load([foldername_fig2{mice},'\neuronIndividuals_new.mat']);
    load([foldername_fig2{mice},'\','circle_results','\','single_cell_firing_profile_S.mat']);
    temp_corr=1-pdist(neuronIndividuals_new{2}.C,'correlation');
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
            rm_corr=corrcoef(cell1,cell2);
            rateMap_corr(ctt)=rm_corr(1,2);
            ctt=ctt+1;
        end
     end
     

%      f1=fit(rateMap_corr',temp_corr','power2');
     subplot(1,3,tk)
     [f1,gof]=fit(rateMap_corr',temp_corr','poly1');
     h1=plot(f1,rateMap_corr(1,1:40:end)',temp_corr(1,1:40:end)');
     set([h1],'color',[0.5 0.5 0.5],'MarkerSize',10)
     gof_collection{tk,1}=gof;
     [rho,p]=corr(rateMap_corr',temp_corr');
     gof_collection{tk,2}=[rho,p]
end

%% new panel g: percentage of place cell
figure;
for tk=1:3  
    mice=mice_all(tk);
    fname=foldername_fig2{mice};
    cluster_file_name=[foldername_fig2{mice},'\','circle_placement_cell_cluster','\','variables_clustering_0.35PC.mat'];
    pc_file_name=[foldername_fig2{mice},'\','circle_info_score_placement_cell','\','place_cells_info_circle_info_scor_binsize10_S_spk.mat'];
  
    neuronName=[foldername_fig2{mice},'\','neuronIndividuals_new.mat'];
    behavName=[foldername_fig2{mice},'\','circle_results\current_condition_behav.mat'];
    threshName=[foldername_fig2{mice},'\','thresh_and_ROI.mat'];
    load(neuronName);
    load(cluster_file_name);
    load(pc_file_name);
    group1=group{2};
    
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
    b1=bar([percentage_pc],'FaceColor','flat');
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

%% 102121 new analysis -- pairwise corr within clust/between clust
% behav
all_behav={};
for i=1:length(foldername_fig2)
    load([foldername_fig2{i},'\','square_Behav.mat']);
    all_behav{i,1}=behav;
    load([foldername_fig2{i},'\','circle_Behav.mat']);
    all_behav{i,2}=behav;
end

% group
gp_rec={};
for i=1:12
    load([foldername_fig2{i},'\','square_placement_cell_cluster','\','variables_clustering_0.35PC.mat'])
    group_model=group{1};

    gp_rec{i,1}=group_model;
    
    load([foldername_fig2{i},'\','circle_placement_cell_cluster','\','variables_clustering_0.35PC.mat'])
    group_model=group{2};

    gp_rec{i,2}=group_model;
end

% for circle
binsize=10;
figure;
for tk=1:12
    mice=mice_all(tk);
    load([foldername_fig2{mice},'\neuronIndividuals_new.mat']);
    [firingrateS,~,~,countTime] = calculatingCellSpatialForSingleData_Suoqin(neuronIndividuals_new{2},all_behav{tk,2}.position,all_behav{tk,2}.time,[0 0 max(all_behav{tk,2}.position,[],1)],binsize,1:size(neuronIndividuals_new{2}.S,1),3*std(neuronIndividuals_new{2}.S,[],2),'S',[],[],[0.1 1000000],10);

    temp_corr=1-pdist(neuronIndividuals_new{2}.C,'correlation');
    dis_dat=pdist(neuronIndividuals_new{2}.centroid);
    rateMap_corr=[];
    
    gp_type=[]; % 1: intra, 0:inter
    ctt=1;
    tic;
    for i=1:length(firingrateS)-1
        for j=i+1:length(firingrateS)
            cell1=firingrateS{i};
            cell2=firingrateS{j};
            cell1=filter2DMatrices(cell1,1);
            cell2=filter2DMatrices(cell2,1);
            fr(countTime==0)=nan; % get rid of non trespassing periods, avoid unnecessary correlation boost

            cell1=reshape(cell1,size(cell1,1)*size(cell1,2),1);
            cell2=reshape(cell2,size(cell2,1)*size(cell2,2),1);
            rm_corr=corrcoef(cell1,cell2,'Rows','complete');
            rateMap_corr(ctt)=rm_corr(1,2);
            
            if(gp_rec{tk,2}(i)==gp_rec{tk,2}(j))
                gp_type(ctt)=1;
            else
                gp_type(ctt)=0;
            end
            ctt=ctt+1;
        end
     end
     toc;

    % plot
    figure;
    for i=1:40:length(rateMap_corr)
     if gp_type(i)==0
         colorr=[126,47,142]/255;
     else
         colorr=[119,172,48]/255;
     end
        
     h1=plot(rateMap_corr(1,i)',temp_corr(1,i)','.','color',colorr); 
     hold on;
    end
    
    figure;
    for i=1:40:length(rateMap_corr)
     if gp_type(i)==1
         colorr=[126,47,142]/255;
     else
         colorr=[119,172,48]/255;
     end
        
     h1=plot(rateMap_corr(1,i)',dis_dat(1,i)','.','color',colorr); 
     hold on;
    end
     
end

