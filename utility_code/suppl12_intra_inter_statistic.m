function [v1,v2,v3,v4,p1,p2,p3]=suppl12_intra_inter_statistic(intra_all_multiGeo,intra_all_multiGeo_tort,intra_shuffle_all_multiGeo_tort,inter_all_multiGeo_tort)

v1=cell2mat(intra_all_multiGeo');
v2=cell2mat(intra_all_multiGeo_tort');
v3=cell2mat(intra_shuffle_all_multiGeo_tort');
v4=cell2mat(inter_all_multiGeo_tort');

v1=distribution_based_resampling(v1,10,10);
v2=distribution_based_resampling(v2,10,10);
v3=distribution_based_resampling(v3,10,10);
v4=distribution_based_resampling(v4,10,10);

disp(['kmeanIntra:', num2str(mean(v1)), '  icaIntra:' ,num2str(mean(v2)), '  icaIntraShuffle:', num2str(mean(v3)), '  icaInter:', num2str(mean(v4))])
disp(['kmeanIntra:', num2str(sem(v1,1)), '  icaIntra:' ,num2str(sem(v2,1)), '  icaIntraShuffle:', num2str(sem(v3,1)), '  icaInter:', num2str(sem(v4,1))])

[h,p1]=kstest2(v1,v2);
[h,p2]=kstest2(v2,v4);
[h,p3]=kstest2(v2,v3);

disp(['kmeanIntra-icaIntra:', num2str(p1), '  icaIntra-icaInter:' ,num2str(p2), '  icaIntra-icaIntraShuffle:', num2str(p3)])

