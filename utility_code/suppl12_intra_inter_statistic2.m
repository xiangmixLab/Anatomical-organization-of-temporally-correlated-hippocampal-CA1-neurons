function [v1,v2,v3,v4,p1,p2,p3]=suppl12_intra_inter_statistic2(intra_all_multiGeo,intra_all_multiGeo_tort,intra_shuffle_all_multiGeo_tort,inter_all_multiGeo_tort)

v1=cell2mat(intra_all_multiGeo');
v2=cell2mat(intra_all_multiGeo_tort');
v3=cell2mat(intra_shuffle_all_multiGeo_tort');
v4=cell2mat(inter_all_multiGeo_tort');

mean(v1)
sem(v1,1)
mean(v2)
sem(v2,1)
mean(v3)
sem(v3,1)
mean(v4)
sem(v4,1)

[h,p1]=kstest2(v2(1:10:end),v4(1:100:end))
[h,p2]=kstest2(v2(1:10:end),v3(1:10000:end))
[h,p3]=kstest2(v1(1:100:end),v2(1:100:end))