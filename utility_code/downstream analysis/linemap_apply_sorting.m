function [fr_sm_111]=linemap_apply_sorting(fr_sm_11,idxx)

for i=1:length(idxx)
    fr_sm_111{i}=fr_sm_11{i}(idxx{i},:);
end

