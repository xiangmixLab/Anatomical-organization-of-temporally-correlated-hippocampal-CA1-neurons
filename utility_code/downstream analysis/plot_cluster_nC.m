function plot_cluster_nC(nC,group)

[~,idx]=sort(group);
imagesc(nC(idx,:));
hold on
        
positionC2 = [0];
for i =1:length(unique(group))
    positionC2(i+1) = positionC2(i)+sum(group == i);
end

for i = 1:length(unique(group))
    plot(1:size(nC,2),ones(1,size(nC,2))*positionC2(i+1),'-','color','r','lineWidth',2);
end

