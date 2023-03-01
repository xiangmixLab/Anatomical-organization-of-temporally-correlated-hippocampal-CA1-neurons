function [firingrateAll_ensem_all,countTime_ensem_all,n_gp_ensem,n_gp]=cluster_ensemble_ratemap_trialPair_analysis_111121(gp1,gp2,clust_to_chk,nnew,behav,binsize)

gp1(gp1~=clust_to_chk)=-1;
gp1(gp1~=-1)=1;

gp2(gp1==-1)=-1;
gp2(gp2~=-1)=gp2(gp2~=-1)-min(gp2(gp2~=-1))+1;

gp_cur={gp1,gp2};

[firingrateAll_ensem_all,countTime_ensem_all,n_gp_ensem,n_gp]=cluster_ensemble_activity_analysis_102821(nnew,gp_cur,behav,binsize);
