function shade_mean_std_plot(data,colorr)
meanVal=mean(data,1);
sem=std(data,[],1);

x = 1:length(meanVal);
x=x-length(meanVal)/2;

y1 = meanVal+sem;
y2 = meanVal-sem;

h = fill([x, fliplr(x)], [y1, fliplr(y2)], colorr);
set(h, 'FaceAlpha', 0.5)
hold on
plot(x,meanVal,'lineWidth',2,'color',colorr);

