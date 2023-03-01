function bdCoor_all=neuronCoor_cal(neuron,varargin)

if isempty(neuron.imageSize)
    neuron.imageSize=size(neuron.Cn);
end

d1=neuron.imageSize(1);
d2=neuron.imageSize(2);

if length(varargin)>0
    thresh=varargin{1};
else
    thresh=0;
end

bdCoor_all={};
for i=1:size(neuron.A,2)
    coor=[];
    
    Ai=full(neuron.A(:,i));

    Ai=Ai*1/max(Ai(:));
    Ait=reshape(Ai,d1,d2);
    Ait=Ait>thresh;
    
    [coor(:,2),coor(:,1)]=find(Ait==1);
    bd=boundary(coor);
    bdCoor=coor(bd,:);
    bdCoor_all{i}=bdCoor;
end


