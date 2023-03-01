function trajectory_firingpos_plot(behav,neuron,behavROI,cell_idx,thresh,markerSizee,temp)

if cell_idx<0
    cell_idx=[1:size(neuron.C,1)];
end

if isequal(temp,'C')
    idx = neuron.C(cell_idx,:)>thresh(cell_idx);
else
    idx = neuron.S(cell_idx,:)>thresh(cell_idx);
end

downsampling = length(neuron.time)/size(neuron.C,2);
if downsampling ~= 1
    neuron.time = double(neuron.time);
    neuron.time = neuron.time(1:downsampling:end);
end
temp = find(diff(behav.time)<=0);
behav.time(temp+1) = behav.time(temp)+1;
neuron.pos = interp1(behav.time,behav.position,neuron.time); %%

plot(neuron.pos(:,1),neuron.pos(:,2),'k')
hold on

if length(cell_idx)>1 % ensemble trace
    idx=mean(idx(:,1:length(neuron.pos(:,1))),1);
    idx=idx>3*std(idx);
end

plot(neuron.pos(idx,1),neuron.pos(idx,2),'r.','MarkerSize',markerSizee)
plot(0,0);
plot(behavROI(1,3),behavROI(1,4));

try % simplified behav do not have obj
    posObjects=round(behav.object);
    if isempty(find(posObjects==0))
        for i5 = 1:size(posObjects,1)
            scatter(posObjects(i5,1),max(max(neuron.pos(:,2)))-posObjects(i5,2)+1,markerSizee,'k','filled')
        end
    end
catch
end

axis ij
shading flat;
axis image
axis off

