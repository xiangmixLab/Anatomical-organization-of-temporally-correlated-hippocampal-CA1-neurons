function peak_rate=suppl18_ensem_individual_ratemap_plots(foldername_fig2,mouse_idx,group_all,behavName,single_field_idx)

peak_rate={};
cond=2;
for tk1=1:length(mouse_idx)
    tk=mouse_idx(tk1);
    load([foldername_fig2{tk},'\','neuronIndividuals_new.mat']);
    load([foldername_fig2{tk},'\',behavName{cond}])
    
    group=group_all{tk,1};       
    uni_group=unique(group);
    uni_group(uni_group==-1)=[];
    
    pct=1;
    figure;
    for i=1:length(uni_group)
        
        % ratemap cal, ensemble and individual
        neuron=neuronIndividuals_new{cond}.copy;
        neuron.delete(setdiff(1:length(size(neuron.C,1)),single_field_idx{tk1,cond}));
%         idx=intersect(find(group==uni_group(i)),single_field_idx{tk1,cond});
        idx=find(group==uni_group(i)); % group should be already processed
        [ensembleMap,firingrateAll,countTime]=cluster_ratemaps(neuron,idx,behav,'S');
               
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
        
        pct=pct+1;
    end
    set(gcf,'renderer','painters');
end
