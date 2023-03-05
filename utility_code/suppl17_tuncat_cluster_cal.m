function [group_LT_tuncat,group_ori_LT_tuncat,group_LT_tuncat_no_rm,group_ori_LT_tuncat_no_rm]=suppl17_tuncat_cluster_cal(dpath,miceNum,trialNum,corr_tun)

all_neurons={};
group_LT_tuncat={};
group_ori_LT_tuncat={};

group_LT_tuncat_no_rm={};
group_ori_LT_tuncat_no_rm={};

t1=tic;
for j=1:miceNum
    load([dpath{j},'\','neuronIndividuals_new_tuncat.mat'])
    all_neurons(j,:)=neuronIndividuals_new_tuncat;
end
disp(['finish loading neuron, ',num2str(toc(t1)),'sec'])

t1=tic;
parfor j=1:miceNum
    for j1=1:trialNum
        n=all_neurons{j,j1}.copy;
        corr_curr=corr_tun{j,j1};
        idx_preserve=corr_curr>0.5; % fur tuncat res, only use neurons with high corr with original trace (probably real neurons)
        n.C=n.C(idx_preserve,:);
        n.S=n.S(idx_preserve,:);

        [g,go]=cluster_determine_by_suoqin_NMF_firstPeakCoph_022422(n,100,10,[]);

        g_full=ones(length(corr_curr),1)*-1;
        g_full(idx_preserve==1)=g;

        go_full=ones(length(corr_curr),1)*-1;
        go_full(idx_preserve==1)=go;
        
        group_LT_tuncat{j,j1}=g_full;
        group_ori_LT_tuncat{j,j1}=go_full;    
        
        [g_no_rm,go_no_rm]=cluster_determine_by_suoqin_NMF_firstPeakCoph_022422(all_neurons{j,j1},100,10,[]);
        group_LT_tuncat_no_rm{j,j1}=g_no_rm;
        group_ori_LT_tuncat_no_rm{j,j1}=go_no_rm;  
        
    end
    disp(['finish, ',num2str(toc(t1)),'sec'])
end