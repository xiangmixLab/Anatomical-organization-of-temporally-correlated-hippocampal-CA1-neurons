function [idxx1]=idx_furtherSelect(idxx,pc)

idxx1=idxx;
for i=1:length(idxx)
    idxx1{i}{1}=intersect(idxx1{i}{1},find(pc==1),'stable');
    idxx1{i}{2}=intersect(idxx1{i}{2},find(pc==1),'stable');
end