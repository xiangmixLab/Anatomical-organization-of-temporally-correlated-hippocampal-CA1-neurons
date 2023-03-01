function vtbt=velo_thresh_by_type_median(all_velo,trials_selection,AD_idx)

all_velo_test=[];
for tk=1:size(all_velo,1)
    for j=trials_selection(tk,:)
        all_velo_test(tk,j)=nanmedian(all_velo{tk,j});
    end
end
all_velo_test(all_velo_test==0)=nan;

all_velo_test_type={};
uni_type=unique(AD_idx);
for i=1:length(uni_type)
    avm=all_velo_test(AD_idx==uni_type(i),:); % Ntg old
    all_velo_test_type{i}=reshape(avm,size(avm,1)*size(avm,2),1);
end

all_velo_test_type_m=fill_nan_to_cellmat(all_velo_test_type);

all_velo_test_type_m{1}(all_velo_test_type_m{1}==0)=nan;

vtbt=nanmean(all_velo_test_type_m{1},1); % Ntg old, AD old, Ntg young, AD young