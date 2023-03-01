% this cluster overlap is still based on rand index, but follows modi et
% al. 2014.
% original rand index = (a+b)/nchoosek(n,2), where a is the number of pairs
% stay together in both and b is the # of pairs separate in both
% This overlap is just (a)/nchoosek(n,2). no separate pairs considered
% 
function [norm_overlap,avg_overlapS,overlap_s]=new_shuf_cluster_overlap(c1,c2)

overlap=new_cluster_overlap(c1,c2);
for i=1:1000
    c1s=c1(randperm(length(c1)));
    c1s=c1s(randperm(length(c1)));
    c2s=c2(randperm(length(c2)));
    c2s=c2s(randperm(length(c2)));
    overlap_s(i)=new_cluster_overlap(c1s,c2s);
end
overlap_s=[overlap_s,overlap];
% norm_overlap=(overlap-min(overlap_s))/(max(overlap_s)-min(overlap_s));
% norm_overlap = prctile(overlap_s,overlap);
% norm_overlap=(overlap-mean(overlap_s))/(std(overlap_s));
% norm_overlap=(overlap-mean(overlap_s))/(max(overlap_s)-mean(overlap_s));
% norm_overlap=overlap/median(overlap_s);
norm_overlap=(overlap-mean(overlap_s))/(1-mean(overlap_s));

avg_overlapS=mean(overlap_s);