function [neuron0,indices_rem]=overlap_handling(neuron)
tic

%% find overlap neurons and delete the one with larger idx
nA=neuron.A;
indices=zeros(size(nA,2),1);
% [KL]=neuron_KL_calculation(neuron);

parfor i=1:size(nA,2)
    dshape=neuron.imageSize;
    A=reshape(nA(:,i),dshape(1),dshape(2))>0;
    sign=overlapDetermine(A,nA,i,dshape,0);
    indices(i)=indices(i)+sign;
end

indices_rem=find(indices>0);
neuron0=neuron.copy;
neuron0.delete(indices_rem)
neuron0.Coor(indices_rem)=[];
toc;