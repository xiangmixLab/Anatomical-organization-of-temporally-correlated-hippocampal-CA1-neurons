function spr=sparsity_calculation(firingrateAllt,cpt)

Po=cpt/sum(cpt(:));
for j=1:length(firingrateAllt)

    if ~isempty(firingrateAllt{j}) && sum(firingrateAllt{j}(:))~=0
        spr(j,1)=sum(sum(firingrateAllt{j}.*Po))^2/sum(sum(firingrateAllt{j}.^2.*Po));
    else
        spr(j,1)=nan;
    end
end