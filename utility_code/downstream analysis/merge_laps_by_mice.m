function fr_bin_11=merge_laps_by_mice(fr_bin_1)
fr_bin_11=cell(1,size(fr_bin_1,2));

for j=1:size(fr_bin_1,2)
    fr_bin_11{1,j}={[],[]};
    for i=1:size(fr_bin_1,1)
        if size(fr_bin_1{i,j}{1},2)~=size(fr_bin_11{1,j}{1},2) && size(fr_bin_11{1,j}{1},2)>0
            fr_bin_1{i,j}{1}=imresize(fr_bin_1{i,j}{1},[size(fr_bin_1{i,j}{1},1),size(fr_bin_11{1,j}{1},2)]);
        end
        if size(fr_bin_1{i,j}{2},2)~=size(fr_bin_11{1,j}{2},2) && size(fr_bin_11{1,j}{2},2)>0
            fr_bin_1{i,j}{2}=imresize(fr_bin_1{i,j}{2},[size(fr_bin_1{i,j}{2},1),size(fr_bin_11{1,j}{2},2)]);
        end
        
        fr_bin_11{1,j}{1}=[fr_bin_11{1,j}{1};fr_bin_1{i,j}{1}];
        fr_bin_11{1,j}{2}=[fr_bin_11{1,j}{2};fr_bin_1{i,j}{2}];       
    end
end
