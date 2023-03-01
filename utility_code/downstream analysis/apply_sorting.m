function [fr_sm_111]=apply_sorting(fr_sm_11,idxx)

for i=1:length(idxx)
    fr_sm_111{i}{1}=fr_sm_11{i}{1}(idxx{i}{1},:);
    fr_sm_111{i}{2}=fr_sm_11{i}{2}(idxx{i}{2},:);
end

