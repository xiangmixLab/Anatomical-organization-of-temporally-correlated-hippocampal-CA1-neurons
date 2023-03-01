function cdf_AD_general(klist)

colorClusters_all=distinguishable_colors(10);

figure;
hold on;
for i=1:length(klist)
    k1=klist{i};

    h1=cdfplot(k1);

    set(h1,'color',colorClusters_all(i,:));

    p1=infer_cdf_loc(k1,nanmedian(k1));

    plot(p1(1),p1(2),'.','color',colorClusters_all(i,:),'MarkerSize',20)
end

