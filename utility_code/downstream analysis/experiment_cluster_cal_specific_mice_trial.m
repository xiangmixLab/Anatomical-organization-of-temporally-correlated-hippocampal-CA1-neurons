function [group,group_ori]=experiment_cluster_cal_specific_mice_trial(dpath,miceNum,trialNum,miceSel,trialSel,K_in)

group=cell(miceNum,trialNum);
group_ori=cell(miceNum,trialNum);
all_neurons=cell(miceNum,trialNum);

%% load neurons
t1=tic;
for j=miceSel
    load([dpath{j},'\','neuronIndividuals_new.mat']);
    for j1=1:size(all_neurons,2)
        all_neurons{j,j1}=Sources2D;
        all_neurons{j,j1}.C=neuronIndividuals_new{j1}.C;
        all_neurons{j,j1}.S=neuronIndividuals_new{j1}.S;
    end
end
disp(['finish loading neuron, ',num2str(toc(t1)),'sec'])

%% generate idx
idxs={};
ctt=1;
for j=miceSel
    for j1=trialSel
        idxs{ctt}=[j,j1];
        ctt=ctt+1;
    end
end

%% parfor cluster calculation
t1=tic;
nums=length(idxs);
g1=cell(nums,1);
g2=cell(nums,1);
parfor j=1:nums
    
    idx_curr=idxs{j};
    
    k=[];
    if ~isempty(K_in)
        if length(K_in)==1
            k=K_in; % single number all mice
        else
            k=K_in(idx_curr(1)); % single number per mice
        end
    end
    
    [g_no_rm,go_no_rm]=cluster_determine_by_suoqin_NMF_firstPeakCoph_022422(all_neurons{idx_curr(1),idx_curr(2)},100,10,k);
    g1{j}=g_no_rm;
    g2{j}=go_no_rm;      
 
    disp(['finish, ',num2str(toc(t1)),'sec'])
end


%% insert clusters
ctt=1;
for j=miceSel
    for j1=trialSel
        idx_curr=idxs{ctt};
        group{idx_curr(1),idx_curr(2)}=g1{ctt};
        group_ori{idx_curr(1),idx_curr(2)}=g2{ctt};
        ctt=ctt+1;
    end
end
