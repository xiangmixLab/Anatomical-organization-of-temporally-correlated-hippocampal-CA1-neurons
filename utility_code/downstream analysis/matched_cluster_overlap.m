function [avg_overlap,preClust_overlap]=matched_cluster_overlap(gp1,gp2)

% suppose gp1, gp2 are two clustering partitions in which the clusts are
% matched, unmatched are labeled with -1
ugp1=unique(gp1);
ugp2=unique(gp2);

ugp1(ugp1==-1)=[];
ugp2(ugp2==-1)=[];

avg_overlap=0;
preClust_overlap=[];

ugp_match=intersect(ugp1,ugp2);
for i=1:length(ugp_match)
    idx1=find(gp1==ugp_match(i));
    idx2=find(gp2==ugp_match(i));
    sameIdx=sum(ismember(idx2,idx1));
    
    preClust_overlap=[preClust_overlap,sameIdx/(length(idx1)+length(idx2)-sameIdx)];
    avg_overlap=avg_overlap+sameIdx/(length(idx1)+length(idx2)-sameIdx);
end

avg_overlap=avg_overlap/length(ugp_match);