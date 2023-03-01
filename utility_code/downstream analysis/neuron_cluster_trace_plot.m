function neuron_cluster_trace_plot(neuron,group2)

colorClusters2=distinguishable_colors(100);

optimalK2=length(unique(group2));

for j = 1:length(unique(group2))
    subplot(length(unique(group2)),1,j);
    rangee=find(group2 == j);
    ctt=1;
    for k=1:length(rangee)
        plot(neuron.C(rangee(k),1:end)'+(ctt-1)*10,'-','color',colorClusters2(j,:));hold on;
        ctt=ctt+1;
    end

end
set(gcf,'renderer','painters');
