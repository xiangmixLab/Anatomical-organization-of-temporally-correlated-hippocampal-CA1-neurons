function shade_line_sem_plot(meanVal,sem,colorr)
x = 1:length(meanVal);
y1 = meanVal+sem;
y2 = meanVal-sem;

plot(meanVal)
h = fill([x, fliplr(x)], [y1, fliplr(y2)], colorr);
set(h, 'FaceAlpha', 0.5)