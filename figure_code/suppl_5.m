%% Supplemental Fig.5: Square experiment, field segregation

% 
run('D:\Xu_clusterting_paper_prep11_2020\final_code\data_prepare\neuron_data_info.m')
all_colors=distinguishable_colors(10);
load(['D:\Xu_clusterting_paper_prep11_2020\final_code\final_cluster_data\cluster_optimal_num\Fig2_cir_rec_clust_original.mat'])

mice_selected=[1,8,11];
condd=1; % square

%% panel A,E,I: cluster footprint
for i=1:length(mice_selected)
    load([foldername_fig2{mice_selected(i)},'\','neuronIndividuals_new.mat'])
    [A_color1,A_color_region1]=DBSCAN_region_quantify_022422(group_ori_fig2{mice_selected(i),condd},neuronIndividuals_new,[]);
    
    subplot(3,2,2*i-1);
    imagesc(A_color1);
    subplot(3,2,2*i);
    imagesc(A_color_region1);
end

%% panel B,F,J: ensemble map
peak_rate={};
for tk=1:3
    mice=mice_selected(tk);
    
    group=group_ori_fig2{mice,condd}; % secondd trial is circle       
    uni_group=unique(group);
    
    load([foldername_fig2{mice_selected(tk)},'\','neuronIndividuals_new.mat'])
    load([foldername_fig2{mice_selected(tk)},'\','behav.mat'])
    
    figure;

    pct=1;

    for i=1:length(unique(group))
        
        % ensemble activity map calculation
        idx=group==uni_group(i);
        [ensembleMap,firingrateAll,countTime]=cluster_ratemaps(neuronIndividuals_new{condd},idx,behavIndividuals{condd},'S');
        t1=double(ensembleMap);
        t1(countTime==0)=nan;
        
        % plot ensemble activity map
        subplot(5,length(unique(group)),pct);
        hold on;
        colormap(jet)
        pcolor(t1);
        axis ij
        shading flat;
        axis off
        peak_rate{tk}(1,i)=max(ensembleMap(:));
        
        % find peak activity map pos
        cen=[];
        for j=1:1:length(firingrateAll)
            tt=firingrateAll{j};
            tt=filter2DMatrices(tt,1);
            [cen(j,1),cen(j,2)]=find(tt==max(tt(:)));
        end
        
        % plot peak positions
        subplot(5,length(unique(group)),pct+length(unique(group)));hold on;
        scatter(round(cen(:,2)),round(cen(:,1)),50,[0.1 0.1 0.1],'.','MarkerFaceAlpha',.5);
        axis ij
        shading flat;
        axis off
       
        
        % find single field cells (the search is random, not 
        [cell_fields,cell_idx]=example_cell_clust_determine(group(idx),i,firingrateAll,countTime,[]);
        
        % plot individual cells
        subplot(5,length(uni_group),pct+length(uni_group)*2);
        ratemap_plot(cell_fields{1},countTime,1,0,[])
        peak_rate{tk}(2,i)=nanmax(filter2DMatrices(cell_fields{1}(:),1));
        axis ij
        shading flat;
        axis off
         
        subplot(5,length(uni_group),pct+length(uni_group)*3);
        ratemap_plot(cell_fields{2},countTime,1,0,[])
        peak_rate{tk}(3,i)=nanmax(filter2DMatrices(cell_fields{2}(:),1));
        axis ij
        shading flat;
        axis off
        
        subplot(5,length(uni_group),pct+length(uni_group)*4);
        ratemap_plot(cell_fields{3},countTime,1,0,[])
        peak_rate{tk}(4,i)=nanmax(filter2DMatrices(cell_fields{3}(:),1));
        axis ij
        shading flat;
        axis off
        
        pct=pct+1;
    end
    set(gcf,'renderer','painters');
end

%% C,G,K: correlation
figure;
gof_collection={};
resample_factor=[40,40,40]; % 

for tk=1:3
    mice=mice_selected(tk);
    
    load([foldername_fig2{mice},'\','neuronIndividuals_new.mat']);
    load([foldername_fig2{mice},'\','behav.mat']);

    % temporal correlation
    temp_corr=1-pdist(neuronIndividuals_new{condd}.C,'correlation');    
    
     % ratemap correlation
    [~,firingrateAll,countTime]=cluster_ratemaps(neuronIndividuals_new{condd},[1:size(neuronIndividuals_new{condd}.C,1)],behavIndividuals{condd},'S'); 
    rateMap_corr=[];
    ctt=1;
    for i=1:length(firingrateAll)-1
        for j=i+1:length(firingrateAll)
            cell1=firingrateAll{i};
            cell2=firingrateAll{j};
            cell1=filter2DMatrices(cell1,1);
            cell2=filter2DMatrices(cell2,1);
            cell1(countTime==0)=0; % get rid of non trespassing periods, avoid unnecessary correlation boost
            cell2(countTime==0)=0;

            cell1=reshape(cell1,size(cell1,1)*size(cell1,2),1);
            cell2=reshape(cell2,size(cell2,1)*size(cell2,2),1);
            rm_corr=corrcoef(cell1,cell2);
            rateMap_corr(ctt)=rm_corr(1,2);
            ctt=ctt+1;
        end
     end
     
     subplot(1,3,tk)
     [f1,gof]=fit(rateMap_corr',temp_corr','poly1');
     h1=plot(f1,rateMap_corr(1,1:resample_factor(tk):end)',temp_corr(1,1:resample_factor(tk):end)');
     ylim([-0.4,1])
     set([h1],'color',[0 0 0],'MarkerSize',10)
     gof_collection{tk,1}=gof;
     [rho,p]=corr(rateMap_corr',temp_corr');
     gof_collection{tk,2}=[rho,p];
     
end
