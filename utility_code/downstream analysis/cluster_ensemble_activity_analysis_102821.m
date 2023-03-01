function [firingrateAll_ensem_all,countTime_ensem_all,n_gp_ensem,n_gp]=cluster_ensemble_activity_analysis_102821(neuron_dat,group_dat,all_behav,binsize)

n_gp={};
n_gp_ensem={};

for tk=1:size(neuron_dat,1)
    for i=1:size(neuron_dat,2)
        n_gp_ensem{tk,i}=Sources2D;
        
        uni_gpdat=unique(group_dat{tk,i});
        uni_gpdat(uni_gpdat==-1)=[];
        
        for j=1:length(uni_gpdat)
            
            ntemp=neuron_dat{tk,i}.copy;
            del_idx=find((group_dat{tk,i}~=uni_gpdat(j))>0);
            ntemp.delete(del_idx);

            n_gp1_C=ntemp.C;
%             n_gp1_S=C_to_peakS(n_gp1_C);
            n_gp1_C_ensem=nanmean(zscore(n_gp1_C,[],2),1);
            n_gp1_C_ensem(n_gp1_C_ensem<0)=0;
            n_gp1_S_ensem=C_to_peakS(n_gp1_C_ensem);

            n_gp_ensem{tk,i}.A=ntemp.A;
            n_gp_ensem{tk,i}.C(j,:)=n_gp1_C_ensem;
            n_gp_ensem{tk,i}.C_raw(j,:)=n_gp1_C_ensem;
            n_gp_ensem{tk,i}.S(j,:)=n_gp1_S_ensem;
            n_gp_ensem{tk,i}.time=ntemp.time;
            
            n_gp{tk,i}{j}=ntemp.copy;
        end
    end
end

% ensemble rate map

firingrateAll_ensem_all={};
countTime_ensem_all={};
for i=1:size(neuron_dat,1)
    for j=1:size(neuron_dat,2)

        [firingrateAll_ensem,~,~,countTime_ensem] = calculatingCellSpatialForSingleData_040321(n_gp_ensem{i,j},all_behav{i,j}.position,all_behav{i,j}.time,[0,0,max(all_behav{i,j}.position,[],1)],binsize,1:size(n_gp_ensem{i,j}.C,1),3*std(n_gp_ensem{i,j}.S,[],2),'S',[],[],[0 inf],10);
     
        firingrateAll_ensem_all{i,j}=firingrateAll_ensem;
        countTime_ensem_all{i,j}=countTime_ensem;
        
    end
end