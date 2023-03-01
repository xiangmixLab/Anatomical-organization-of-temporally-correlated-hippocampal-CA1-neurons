function [group_shared,group12,group21,group12t,group21t]=determineSharedCells_new_070121(group1,group2)
uni_g1=unique(group1);
uni_g2=unique(group2);

group12=group1*0; % use group1's idx to replace group2's element
group21=group1*0; % use group2's idx to replace group1's element

group12t=group1*0;
group21t=group1*0;

g1_clust_in_g2=[];
uv_all={};
% for i=1:length(uni_g1)
%     idxx=find(group1==uni_g1(i));
%     g2_idxx=group2(idxx);
%     g1_clust_in_g2(i)=[mode(g2_idxx)];
% end

for i=1:length(uni_g1)
    idxx=find(group1==uni_g1(i));
    g2_idxx=group2(idxx);
    
    uv=histc(g2_idxx, unique(g2_idxx));
    uv_mat=zeros(1,max(group1));
    uv_mat(unique(g2_idxx))=uv;
    uv_all{i,1}=(uv_mat./length(g2_idxx)); % ratio
    g1_clust_in_g2(i)=min(find(uv_all{i,1}==max(uv_all{i,1})));
end

while true
    if length(g1_clust_in_g2)==length(unique(g1_clust_in_g2)) % looping until all group has been reassigned
        break;
    else % if there's duplicated assignment, compare the composition ratio and choose the largest one
        for i=1:length(uni_g1)
            idx_same=find(g1_clust_in_g2==i); %
            if length(idx_same)>1
                nnum_compare=cell2mat(uv_all(idx_same));
                nnum_compare1=nnum_compare(:,i);
                max_num=max(nnum_compare1,[],1);
                max_choice=find(max_num==max(max_num));
                if length(max_choice)==1
                    idx_same(max_choice)=[];
                    for l=1:length(idx_same)
                        uv_all{idx_same(l)}(i)=-1;
                        g1_clust_in_g2(idx_same(l))=min(find(uv_all{idx_same(l)}==max(uv_all{idx_same(l)})));
                    end
                else
                    max_choice=max(max_choice); % if two g1 clust have same largest composition ratio in a g2 clust, the clust will be assigned to the relatively larger idx
                    idx_same(max_choice)=[];
                    for l=1:length(idx_same)
                        uv_all{idx_same(l)}(i)=-1;
                        g1_clust_in_g2(idx_same(l))=min(find(uv_all{idx_same(l)}==max(uv_all{idx_same(l)})));
                    end
                end
            end
        end
    end
end

    

for i=1:length(uni_g1)
    group12t(logical((group2==g1_clust_in_g2(i)).*(group1==uni_g1(i))))=i;  
    group12(logical((group2==g1_clust_in_g2(i))))=i;
end


group_shared=group12t>0;