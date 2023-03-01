function [clustNum_perMice,clustNum_allMice]=cluster_num_aggregation(group)

clustNum=[];
ctt=1;
for l=1:size(group,1)
    for j1=1:size(group,2)
        clustNum(l,j1)=max(group{l,j1});
        ctt=ctt+1;
    end
end
clustNum_perMice=floor(mean(clustNum,2));
clustNum_allMice=floor(mean(clustNum(:)));