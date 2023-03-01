function shade_mean_sem_plot(data,colorr)
meanVal=mean(data,1);
sem=std(data,[],1)/size(data,1).^0.5;

x = 1:length(meanVal);
y1 = meanVal+sem;
y2 = meanVal-sem;

plot(meanVal,'lineWidth',2,'color',colorr);
h = fill([x, fliplr(x)], [y1, fliplr(y2)], colorr);
set(h, 'FaceAlpha', 0.5)