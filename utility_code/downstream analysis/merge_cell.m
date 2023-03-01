function dat_merge=merge_cell(dat)

dat_merge=[];
for i=1:length(dat)
    dat_merge=[dat_merge;dat{i}];
end