function [Asize1,ind_del,neuron0]=neuronAreaAndTrim(neuron,threshh)
Asize1=[];
for j=1:size(neuron.A,2)
    Atemp=reshape(neuron.A(:,j),neuron.imageSize(1),neuron.imageSize(2));
    Atemp=Atemp>0;
    Asize1(j,1)=sum(Atemp(:));
end

ind_del=find(Asize1>=threshh);
neuron0=neuron.copy;
neuron0.delete(ind_del);
neuron0.Coor=neuronCoor_cal(neuron0,0.5);