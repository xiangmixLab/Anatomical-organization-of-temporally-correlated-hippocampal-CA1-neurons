%% Supplemental Fig.17: targeted NMF

%% we use linear track and Ai163 datasets to try out

%% 0. file pathes (neuron videos)
multigeo_data_path={
    'D:\Xu_clusterting_paper_prep11_2020\final data\final_neuron_behav_data\Fig_1_3_tri_cir_sqr\M3411';		
    'D:\Xu_clusterting_paper_prep11_2020\final data\final_neuron_behav_data\Fig_1_3_tri_cir_sqr\M3412';		
    'D:\Xu_clusterting_paper_prep11_2020\final data\final_neuron_behav_data\Fig_1_3_tri_cir_sqr\M3421F';		
    'D:\Xu_clusterting_paper_prep11_2020\final data\final_neuron_behav_data\Fig_1_3_tri_cir_sqr\M3422F';		
    'D:\Xu_clusterting_paper_prep11_2020\final data\final_neuron_behav_data\Fig_1_3_tri_cir_sqr\M3424F';		
    'D:\Xu_clusterting_paper_prep11_2020\final data\final_neuron_behav_data\Fig_1_3_tri_cir_sqr\M3425F';		
}

cir_rec_data_path={
    'D:\Xu_clusterting_paper_prep11_2020\final data\final_neuron_behav_data\Fig_2_cir_rec\M3411';	
    'D:\Xu_clusterting_paper_prep11_2020\final data\final_neuron_behav_data\Fig_2_cir_rec\M3412';
    'D:\Xu_clusterting_paper_prep11_2020\final data\final_neuron_behav_data\Fig_2_cir_rec\M3413';		
    'D:\Xu_clusterting_paper_prep11_2020\final data\final_neuron_behav_data\Fig_2_cir_rec\M3414';
    'D:\Xu_clusterting_paper_prep11_2020\final data\final_neuron_behav_data\Fig_2_cir_rec\M3415';		
    'D:\Xu_clusterting_paper_prep11_2020\final data\final_neuron_behav_data\Fig_2_cir_rec\M3421F';
    'D:\Xu_clusterting_paper_prep11_2020\final data\final_neuron_behav_data\Fig_2_cir_rec\M3422F';		
    'D:\Xu_clusterting_paper_prep11_2020\final data\final_neuron_behav_data\Fig_2_cir_rec\M3423F';	
    'D:\Xu_clusterting_paper_prep11_2020\final data\final_neuron_behav_data\Fig_2_cir_rec\M3424F';		
    'D:\Xu_clusterting_paper_prep11_2020\final data\final_neuron_behav_data\Fig_2_cir_rec\M3425F';	
    'D:\Xu_clusterting_paper_prep11_2020\final data\final_neuron_behav_data\Fig_2_cir_rec\M3426F';		
    'D:\Xu_clusterting_paper_prep11_2020\final data\final_neuron_behav_data\Fig_2_cir_rec\M3427F';			
}

linearTrack_data_path={
    'D:\Xu_clusterting_paper_prep11_2020\final data\final_neuron_behav_data\Fig_3_linear_track\M3411'		
    'D:\Xu_clusterting_paper_prep11_2020\final data\final_neuron_behav_data\Fig_3_linear_track\M3412'		
    'D:\Xu_clusterting_paper_prep11_2020\final data\final_neuron_behav_data\Fig_3_linear_track\M3421F'		
    'D:\Xu_clusterting_paper_prep11_2020\final data\final_neuron_behav_data\Fig_3_linear_track\M3422F'		
    'D:\Xu_clusterting_paper_prep11_2020\final data\final_neuron_behav_data\Fig_3_linear_track\M3424F'		
    'D:\Xu_clusterting_paper_prep11_2020\final data\final_neuron_behav_data\Fig_3_linear_track\M3425F'		
}
barrier_data_path={
    'D:\Xu_clusterting_paper_prep11_2020\final data\final_neuron_behav_data\Fig_3_barrier\101F'	
    'D:\Xu_clusterting_paper_prep11_2020\final data\final_neuron_behav_data\Fig_3_barrier\102F'	
    'D:\Xu_clusterting_paper_prep11_2020\final data\final_neuron_behav_data\Fig_3_barrier\103F'		
    'D:\Xu_clusterting_paper_prep11_2020\final data\final_neuron_behav_data\Fig_3_barrier\2833M'	
    'D:\Xu_clusterting_paper_prep11_2020\final data\final_neuron_behav_data\Fig_3_barrier\840F'		
}

%% 1. generate .h5 files
suppl17_h5_gen(linearTrack_data_path);
suppl17_h5_gen(barrier_data_path);
suppl17_h5_gen(multigeo_data_path);
suppl17_h5_gen(cir_rec_data_path);

%% 2. generate mask
suppl17_mask_gen(linearTrack_data_path)
suppl17_mask_gen(barrier_data_path)
suppl17_mask_gen(multigeo_data_path)
suppl17_mask_gen(cir_rec_data_path)

%% 3. identify false cells, bg
% the paper indicate for one-photon vid, the algorithm use cleaned CNMF-E
% mask and generated unmixed traces for each detected area

% handled in TUnCaT python code

%% 4. generate new neuronIndividuals_new
suppl17_CNMFE_style_data_generation(linearTrack_data_path)
suppl17_CNMFE_style_data_generation(barrier_data_path)
suppl17_CNMFE_style_data_generation(multigeo_data_path)
suppl17_CNMFE_style_data_generation(cir_rec_data_path)

%% 4.1 check TUnCaT quality: exclude badly detected neurons 
preserve_idxs=exclude_neurons_by_temporal(linearTrack_data_path);

%% 4.2 use preserved neuron idx to generate new CNMF-E and TUnCaT dataset
suppl17_trimmed_neuron(linearTrack_data_path,preserve_idxs)

%% 5.CNMF-E and TUnCaT trace correlation
[corr_all_LT,corr_per_mice_LT]=suppl17_tuncat_CNMFE_corr(linearTrack_data_path);

%% 6. CNMF-E and TUnCaT trace above-thres spike diff, perprotion
[peakDiff_LT,peakDiff_LT_percentage]=suppl17_tuncat_CNMFE_peakTotalNumDiff(linearTrack_data_path);
[peakDiff_LT,peakDiff_LT_percentage]=suppl17_tuncat_CNMFE_peakDiff(linearTrack_data_path); % figure generation

%% 7.  calculate tuncat clusters
[group_LT_tuncat_no_rm,group_ori_LT_tuncat_no_rm]=experiment_cluster_cal_tuncat(linearTrack_data_path,1,1,[]);

[group_barrier_tuncat,group_ori_barrier_tuncat,group_barrier_tuncat_no_rm,group_ori_barrier_tuncat_no_rm]=suppl17_tuncat_cluster_cal(barrier_data_path,5,6);
[group_mgeo_tuncat,group_ori_mgeo_tuncat,group_mgeo_tuncat_no_rm,group_ori_mgeo_tuncat_no_rm]=suppl17_tuncat_cluster_cal(multigeo_data_path,6,6);
[group_cr_tuncat,group_ori_cr_tuncat,group_cr_tuncat_no_rm,group_ori_cr_tuncat_no_rm]=suppl17_tuncat_cluster_cal(cir_rec_data_path,12,2);

%% 7.1  calculate tuncat cluster, different num, with original, size
[group_LT_tuncat_no_rm_k,group_ori_LT_tuncat_no_rm_k]=suppl17_tuncat_cluster_cal_2(linearTrack_data_path,6,1);
[group_LT_k,group_ori_LT_k]=suppl17_cluster_cal_3(linearTrack_data_path,6,1);

[A_color_LT_tun,A_color_region_LT_tun,avg_region_LT_tun,avg_region_shuf_LT_tun]=suppl17_cluster_region_cal_tuncat(group_ori_LT_tuncat_no_rm_k,linearTrack_data_path);
[A_color_LT,A_color_region_LT,avg_region_LT,avg_region_shuf_LT]=suppl17_cluster_region_cal(group_ori_LT_k,linearTrack_data_path);


%% 8. illustration
%%%%%%%%%%%% original neuron 
neuronIdx={[4,5,7],[4:6],[4:6],[4:6],[4:6],[4,8,6]}
figure;
for j=1:6
    subplot(1,6,j);
    load([linearTrack_data_path{j},'\','neuronIndividuals_new.mat'])
    imagesc(neuronIndividuals_new{1}.Cn);hold on;
    for j1=1:length(neuronIdx{j})
        plot(neuronIndividuals_new{1}.Coor{neuronIdx{j}(j1)}(1,2:end),neuronIndividuals_new{1}.Coor{neuronIdx{j}(j1)}(2,2:end),'-','color','r','lineWidth',2);
        text(mean(neuronIndividuals_new{1}.Coor{neuronIdx{j}(j1)}(1,2:end)),mean(neuronIndividuals_new{1}.Coor{neuronIdx{j}(j1)}(2,2:end)),num2str(j1));
    end
end

%%%%%%%%%%%% example traces to show bg sub, WITH THRESHOLDS
figure;hold on
for j=1
%     subplot(1,6,j);hold on;
    load([linearTrack_data_path{j},'\','neuronIndividuals_new.mat'])
    factor=100;
    nC=zscore(neuronIndividuals_new{1}.C,[],2);
    Cpeaks=C_to_peakS(nC);
    
    for j1=1:length(neuronIdx{j})
        plot(nC(neuronIdx{j}(j1),:)+factor,'color','b');
        thres=3*std(Cpeaks(neuronIdx{j}(j1),:),[],2);
        line([0,size(nC,2)],[thres,thres]+factor,'lineStyle','--','color','r');
        factor=factor-30;
    end
    set(gcf,'renderer','painters');
end

figure;hold on
for j=1
%     subplot(1,6,j);hold on;
    load([linearTrack_data_path{j},'\','neuronIndividuals_new_tuncat.mat'])
    factor=100;
    
    nC=zscore(neuronIndividuals_new_tuncat{1}.C,[],2);
    Cpeaks=C_to_peakS(nC);
    
    for j1=1:length(neuronIdx{j})
        plot(nC(neuronIdx{j}(j1),:)+factor,'color','g')
        thres=3*std(Cpeaks(neuronIdx{j}(j1),:),[],2);
        line([0,size(nC,2)],[thres,thres]+factor,'lineStyle','--','color','r');

        factor=factor-30;
    end
    set(gcf,'renderer','painters');
end

%%%%%%%%%%%% TUN footprint
subplot(121)
imagesc(A_color_LT_tun{1}{max(group_ori_LT_tuncat_no_rm{1})})
subplot(122)
imagesc(A_color_region_LT_tun{1}{max(group_ori_LT_tuncat_no_rm{1})})

%%%%%%%%%%% original footprint
subplot(121)
imagesc(A_color_LT{1}{max(group_ori_LT{1})})
subplot(122)
imagesc(A_color_region_LT{1}{max(group_ori_LT{1})})

%%%%%%%%%%%% TUN cluster size
cond=1;
avg_region_cond1_mat=[];
for i=1:size(avg_region_LT_tun,1)
    avg_region_cond1_mat(i,:)=cell2mat(avg_region_LT_tun{i,cond}); 
end

avg_region_cond1_shuf_mat=[];
for i=1:size(avg_region_shuf_LT_tun,1) 
    shuf_region=avg_region_shuf_LT_tun{i,cond};
    shuf_region_mat=[];
    for j=1:length(shuf_region)
        shuf_region_mat(j)=mean(shuf_region{j});
    end
    avg_region_cond1_shuf_mat(i,:)=shuf_region_mat;
end

mean_avg_max_region_mat=mean(avg_region_cond1_mat,1);
se_avg_max_region_mat=std(avg_region_cond1_mat,[],1)/(size(avg_region_cond1_mat,1)^0.5);

mean_avg_max_reg_shuf_all_mat=mean(avg_region_cond1_shuf_mat(:,2:end),1);
se_avg_max_reg_shuf_all_mat=std(avg_region_cond1_shuf_mat(:,2:end),[],1)/(size(avg_region_cond1_shuf_mat(:,2:end),1)^0.5);

errorbar(0:length(mean_avg_max_region_mat)-1,mean_avg_max_region_mat,se_avg_max_region_mat,'-o','color',[0    0.4470    0.7410]);
hold on;
errorbar(0:length(mean_avg_max_reg_shuf_all_mat)-1,mean_avg_max_reg_shuf_all_mat,se_avg_max_reg_shuf_all_mat,'-s','color',[0.8500    0.3250    0.0980]);

signrank(mean_avg_max_region_mat,mean_avg_max_reg_shuf_all_mat)

p=[];
for i=1:size(avg_region_cond1_mat,2)
    [p(i,1)]=signrank(avg_region_cond1_mat(:,i),avg_region_cond1_shuf_mat(:,i+1));
end
p
%%%%%%%%%%%% original cluster size
cond=1;
avg_region_cond1_mat=[];
for i=1:size(avg_region_LT_tun,1)
    avg_region_cond1_mat(i,:)=cell2mat(avg_region_LT{i,cond}); 
end

avg_region_cond1_shuf_mat=[];
for i=1:size(avg_region_shuf_LT,1) 
    shuf_region=avg_region_shuf_LT{i,cond};
    shuf_region_mat=[];
    for j=1:length(shuf_region)
        shuf_region_mat(j)=mean(shuf_region{j});
    end
    avg_region_cond1_shuf_mat(i,:)=shuf_region_mat;
end

mean_avg_max_region_mat=mean(avg_region_cond1_mat,1);
se_avg_max_region_mat=std(avg_region_cond1_mat,[],1)/(size(avg_region_cond1_mat,1)^0.5);

mean_avg_max_reg_shuf_all_mat=mean(avg_region_cond1_shuf_mat(:,2:end),1);
se_avg_max_reg_shuf_all_mat=std(avg_region_cond1_shuf_mat(:,2:end),[],1)/(size(avg_region_cond1_shuf_mat(:,2:end),1)^0.5);

errorbar(0:length(mean_avg_max_region_mat)-1,mean_avg_max_region_mat,se_avg_max_region_mat,'-o','color',[0    0.4470    0.7410]);
hold on;
errorbar(0:length(mean_avg_max_reg_shuf_all_mat)-1,mean_avg_max_reg_shuf_all_mat,se_avg_max_reg_shuf_all_mat,'-s','color',[0.8500    0.3250    0.0980]);

signrank(mean_avg_max_region_mat,mean_avg_max_reg_shuf_all_mat)

p=[];
for i=1:size(avg_region_cond1_mat,2)
    [p(i,1)]=signrank(avg_region_cond1_mat(:,i),avg_region_cond1_shuf_mat(:,i+1));
end
p