function [neuron]=thresholdNeuronFootprint(neuron,thresh)

if isempty(neuron.imageSize)
    neuron.imageSize=size(neuron.Cn);
end

d1=neuron.imageSize(1);
d2=neuron.imageSize(2);

if isempty(thresh)
    thresh=0;
end

neuron1=neuron.copy;

Ai=full(neuron.A);
Ai=Ai*1/max(Ai(:));
Ai(Ai<thresh)=0;
neuron1.A=Ai;

