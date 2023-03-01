function plot_contour_Proj(vidName,neuronName)

[~, ~, fExt] = fileparts(vidName);
if strcmp(lower(fExt),'.avi')
    V=general_avi_read(vidName);
    V=cat(3,V{:});
else
    if strcmp(lower(fExt),'.mat')
        V=importdata(vidName);
    end
end

load(neuronName);
bdCoor_all=neuronCoor_cal(neuron,0.5);

imagesc(neuron.Cn); hold on;
for i=1:length(neuron.Coor)
    plot(bdCoor_all{i}(:,1),bdCoor_all{i}(:,2),'r','lineWidth',2);
end