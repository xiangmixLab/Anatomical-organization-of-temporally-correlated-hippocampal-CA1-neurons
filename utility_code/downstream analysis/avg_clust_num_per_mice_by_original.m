function clustNum=avg_clust_num_per_mice_by_original(gp_ori)

clustNum=[];
for l=1:size(gp_ori,1)
    for l1=1:size(gp_ori,2)
        clustNum(l,l1)=max(gp_ori{l,l1});
    end
end
clustNum=round(mean(clustNum,2));