% suppl 12-13: place cells, 
run('D:\Xu_clusterting_paper_prep11_2020\final_code\data_prepare\neuron_data_info.m')

%% 1. place cell calculation
foldername=foldername_fig2;

all_behav={};
for i=1:length(foldername)
    load([foldername{i},'\','behav.mat']);
    for j=1:size(behavIndividuals,2) 
        all_behav{i,j}=behavIndividuals{j};
        all_behav{i,j}.VidObj=[];
    end
end

% PC
all_pc=cell(12,2);
all_infoscore=cell(12,2);
all_infoscore_norm=cell(12,2);
all_coherence=cell(12,2);
tic;
for i=1:12
    load([foldername{i},'\','neuronIndividuals_new.mat']);
    for j=1:2
    
        [place_cells,infoScore,infoScore_norm,coherencee] = permutingSpike_adapt_040821(neuronIndividuals_new{j},all_behav{i,j}.position,all_behav{i,j}.time,'S',0,10,5,'all',0.4);  

        all_pc{i,j}=place_cells;
        all_infoscore{i,j}=infoScore;
        all_infoscore_norm{i,j}=infoScore_norm;
        all_coherence{i,j}=coherencee;

    end
end
toc;

%% config 
cond=2; % suppl 12 use 2 (circle), suppl 13 use 1 (square)
resample_factor=[2,1,1]; % visualization dot resampling

mouse_idx=[1,8,11];
color_clust=distinguishable_colors(10);
%% anatomical footprint, place cell
figure;
for tk1=1:length(mouse_idx)
    tk=mouse_idx(tk1);
    load([foldername_fig2{tk},'\','neuronIndividuals_new.mat']);
    
    gp=group_ori_fig2{tk,cond};
    pc_idx=all_pc{tk,cond}{2};
    gp(setdiff([1:length(gp)],pc_idx))=-1;
    A_color=cluster_spatial_footprint_colormap(neuronIndividuals_new,d1,d2,color_clust,gp,0.75);

    subplot(1,length(mouse_idx),tk1)
    imagesc(A_color);
end
%% ensemble maps, place cell
peak_rate={};
single_cell_idx={};
for tk1=1:length(mouse_idx)
    tk=mouse_idx(tk1);
    load([foldername_fig2{tk},'\','neuronIndividuals_new.mat']);
    
    group=group_ori_fig2{tk,cond};       
    uni_group=unique(group);
    uni_group(uni_group==-1)=[];
    
    pct=1;
    figure;
    for i=1:length(uni_group)
        
        % ratemap cal, ensemble and individual
        neuron=neuronIndividuals_new{cond}.copy;
        idx=intersect(find(group==uni_group(i)),all_pc{tk,cond}{2});
%         idx=find(group==uni_group(i));
        [ensembleMap,firingrateAll,countTime]=cluster_ratemaps(neuron,idx,all_behav{tk,cond},'S');
               
        % find individual peak position
        max_inds={};
        prates={};
        PF_radius={};
        for j=1:1:length(firingrateAll)
            tt=firingrateAll{j};
            tt(countTime==0)=nan;
            tt=filter2DMatrices(tt,1);
            [max_inds{j}(1),max_inds{j}(2)]=find(tt==max(tt(:)));  % find all local peak locations, good for multi-field cells
        end

        % ensemble ratemap plot
        subplot(4,length(uni_group),pct);
        ratemap_plot(ensembleMap,countTime,1,0,[])
        hold on;
        colormap(jet)
        
        for j=1:1:length(idx) % plot individual peaks
            cen=max_inds{j}; 
            scatter1=scatter(round(cen(1,2)),round(cen(1,1)),50,[0.1 0.1 0.1],'.','MarkerFaceAlpha',.5);
        end
        
        % record peak rate for image illustration
        peak_rate{tk}(1,i)=nanmax(ensembleMap(:));
        
        % find single field cells
        [cell_fields,cell_idx]=single_field_cell_clust_determine(group(idx),i,firingrateAll,countTime,[]);
        
        % plot individual cells
        subplot(4,length(uni_group),pct+length(uni_group));
        ratemap_plot(cell_fields{1},countTime,1,0,[])
        peak_rate{tk}(2,i)=nanmax(cell_fields{1}(:));
         
        subplot(4,length(uni_group),pct+length(uni_group)*2);
        ratemap_plot(cell_fields{2},countTime,1,0,[])
        peak_rate{tk}(3,i)=nanmax(cell_fields{2}(:));
        
        subplot(4,length(uni_group),pct+length(uni_group)*3);
        ratemap_plot(cell_fields{3},countTime,1,0,[])
        peak_rate{tk}(4,i)=nanmax(cell_fields{3}(:));
        
        single_cell_idx{tk}{i}=cell_idx;
        pct=pct+1;
    end
    set(gcf,'renderer','painters');
end

%% traj-firing dot maps, place cell
for tk1=1:length(mouse_idx)
    tk=mouse_idx(tk1);
    load([foldername_fig2{tk},'\','neuronIndividuals_new.mat']);
    
    group=group_ori_fig2{tk,cond};       
    uni_group=unique(group);
    uni_group(uni_group==-1)=[];
    
    pct=1;
    figure;
    for i=1:length(uni_group)
        
        % ratemap cal, ensemble and individual
        neuron=neuronIndividuals_new{cond}.copy;
        idx=intersect(find(group==uni_group(i)),all_pc{tk,cond}{2});
        
        behav=all_behav{tk,cond};
        behavROI=all_behav{tk,cond}.ROI;
        
        n1=Sources2D;
        n1.C=neuron.C(idx,:);
        n1.S=neuron.S(idx,:);
        n1.time=neuron.time;
        
        ne=Sources2D;
        ne.C=mean(n1.C,1);
        ne.S=mean(n1.S,1);
        ne.time=n1.time;
        
        % ensemble ratemap plot
        subplot(4,length(uni_group),pct);    
        trajectory_firingpos_plot(behav,ne,behavROI,1,3*std(ne.S,[],2),10,'S');
        
        cellIdx_cur=single_cell_idx{tk}{i};
        
        % plot individual cells
        subplot(4,length(uni_group),pct+length(uni_group));
        trajectory_firingpos_plot(behav,n1,behavROI,cellIdx_cur(1),3*std(n1.S,[],2),10,'S');
         
        subplot(4,length(uni_group),pct+length(uni_group)*2);
        trajectory_firingpos_plot(behav,n1,behavROI,cellIdx_cur(2),3*std(n1.S,[],2),10,'S');
        
        subplot(4,length(uni_group),pct+length(uni_group)*3);
        trajectory_firingpos_plot(behav,n1,behavROI,cellIdx_cur(3),3*std(n1.S,[],2),10,'S');
        
        pct=pct+1;
    end
    set(gcf,'renderer','painters');
end

%% correlation plot
figure;
gof_collection={};
for tk=1:3
    mice=mouse_idx(tk);
    
    load([foldername_fig2{tk},'\','neuronIndividuals_new.mat']);
    
    neuron=neuronIndividuals_new{cond}.copy;
    idx=all_pc{tk,cond}{2};
        
    group=group_ori_fig2{tk,cond};       
    uni_group=unique(group);
    uni_group(uni_group==-1)=[];
    
    temp_corr=1-pdist(neuronIndividuals_new{cond}.C(idx,:),'correlation');
    
    [~,firingrateAll,countTime]=cluster_ratemaps(neuron,idx,all_behav{tk,cond},'S');
    
    rateMap_corr=[];
    ctt=1;
    
    
    for i=1:length(firingrateAll)-1
        for j=i+1:length(firingrateAll)
            cell1=firingrateAll{i};
            cell2=firingrateAll{j};
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
     
     subplot(1,3,tk)
     [f1,gof]=fit(rateMap_corr',temp_corr','poly1');
     h1=plot(f1,rateMap_corr(1,1:resample_factor(tk):end)',temp_corr(1,1:resample_factor(tk):end)');
     ylim([-0.4,1])
     set([h1],'color',[0 0 0],'MarkerSize',10)
     gof_collection{tk,1}=gof;
     [rho,p]=corr(rateMap_corr(1,1:resample_factor(tk):end)',temp_corr(1,1:resample_factor(tk):end)');
     gof_collection{tk,2}=[rho,p]
     
end

