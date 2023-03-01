function [idxx1]=linemap_idx_furtherSelect(idxx,pc)

idxx1=idxx;
for i=1:length(idxx)
    idxx1{i}=intersect(idxx1{i},find(pc==1),'stable');
end