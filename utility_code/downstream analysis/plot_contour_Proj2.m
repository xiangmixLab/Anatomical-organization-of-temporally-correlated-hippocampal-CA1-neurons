function plot_contour_Proj2(neuron)


bdCoor_all=neuronCoor_cal(neuron,0.5);

imagesc(neuron.Cn); hold on;
for i=1:length(bdCoor_all)
    plot(bdCoor_all{i}(:,1),bdCoor_all{i}(:,2),'r','lineWidth',2);
end