function displayHeatmap_general_031721(nC,group2)
perm=[];
for j = 1:length(unique(group2))
    perm = [perm;find(group2 == j)];
end
CM1=squareform(1-pdist(nC,'correlation'));
CM1(logical(eye(size(CM1)))) = 1;
if 0
    imagesc(flip(CM1(perm,perm)));
else
    imagesc(CM1(perm,perm));
end
