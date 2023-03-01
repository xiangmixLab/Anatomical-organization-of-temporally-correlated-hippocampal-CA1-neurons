function plot_cluster_nC_trace(nC,group)

[~,idx]=sort(group);
% 
figure;
uni_gp=unique(group);
colorss=distinguishable_colors(10);
for i=1:length(uni_gp)
    subplot(length(uni_gp),1,i);
    nC_gp=nC(group==uni_gp(i),:);
    for k=1:2:size(nC_gp,1)
        plot(zscore(nC_gp(k,:))+(k-1)*2,'color',colorss(i,:));hold on;
    end
    axis off
end
