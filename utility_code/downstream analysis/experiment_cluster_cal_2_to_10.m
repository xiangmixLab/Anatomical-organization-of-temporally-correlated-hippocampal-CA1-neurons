function [group,group_ori]=experiment_cluster_cal_2_to_10(dpath,miceNum,trialNum)

group=cell(miceNum,trialNum);
group_ori=cell(miceNum,trialNum);
all_neurons=cell(miceNum,trialNum);

%% load neurons
t1=tic;
for j=1:miceNum
    load([dpath{j},'\','neuronIndividuals_new.mat']);
    for j1=1:size(all_neurons,2)
        all_neurons{j,j1}=Sources2D;
        all_neurons{j,j1}.C=neuronIndividuals_new{j1}.C;
        all_neurons{j,j1}.S=neuronIndividuals_new{j1}.S;
        group{j,j1}={};
        group_ori{j,j1}={};
    end
end
disp(['finish loading neuron, ',num2str(toc(t1)),'sec'])

%% generate idx
idxs={};
ctt=1;
for j=1:miceNum
    for j1=1:trialNum
        for k=2:10
            idxs{ctt}=[j,j1,k];
            ctt=ctt+1;
        end
    end
end

%% parfor cluster calculation
t1=tic;
nums=length(idxs);
g1=cell(nums,1);
g2=cell(nums,1);
parfor j=1:nums
    
    idx_curr=idxs{j};
    [g_no_rm,go_no_rm]=cluster_determine_by_suoqin_NMF_firstPeakCoph_022422(all_neurons{idx_curr(1),idx_curr(2)},100,10,idx_curr(3));
    g1{j}=g_no_rm;
    g2{j}=go_no_rm;      
 
    disp(['finish, ',num2str(toc(t1)),'sec'])
end


%% insert clusters
ctt=1;
for j=1:miceNum
    for j1=1:trialNum
        for k=2:10
            idx_curr=idxs{ctt};
            group{idx_curr(1),idx_curr(2)}{idx_curr(3)}=g1{ctt};
            group_ori{idx_curr(1),idx_curr(2)}{idx_curr(3)}=g2{ctt};
            ctt=ctt+1;
        end
    end
end
