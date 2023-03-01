function cluster_colorbar(group)
color_bar=[];
color_all=distinguishable_colors(10);
perm=[];
for j = 1:length(unique(group))
    perm = [perm;find(group == j)];
end
for i=1:length(perm)
    color_bar(i,1,:)=color_all(group(perm(i)),:);
end
imagesc(color_bar);
