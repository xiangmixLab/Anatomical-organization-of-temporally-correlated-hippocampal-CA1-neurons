%% simple graph cell register by pairwise dis

c1=neurons{1}.centroid;
c2=neurons{2}.centroid;

aa=squareform(pdist(c1));
bb=squareform(pdist(c2));

M=[];

for i=1:size(aa,1)
    for j=1:size(bb,1)
        neighboor_c1_i=aa(i,:);
        neighboor_c2_i=bb(j,:);
        
        neighboor_c1_i=sort(neighboor_c1_i);
        neighboor_c2_i=sort(neighboor_c2_i);
        
        size_choose=round(min(length(neighboor_c1_i)*0.4,length(neighboor_c2_i)*0.4));
        
        neighboor_c1_i=neighboor_c1_i(1:size_choose);
        neighboor_c2_i=neighboor_c2_i(1:size_choose);
        
%         corr_t=corrcoef(neighboor_c1_i,neighboor_c2_i);
%         M(i,j)=corr_t(2);
        M(i,j)= nanmean((neighboor_c2_i-neighboor_c1_i).^2);
    end
end

M_99tile=quantile(M(:),0.5);

cell_register=[];
for i=1:size(M,1)
    cell_register(i,1)=i;
%     [~,idx]=sort(M(i,:),'descend'); % large to small
    [~,idx]=sort(M(i,:));

    for j=1:length(idx)
        max_deviation=15;
        c1t=c1(i,:);
        c2t=c2(idx(j),:);
        distt=sum((c1t-c2t).^2)^0.5;
        if distt<max_deviation
            cell_register(i,2)=idx(j);
            break;
        end
    end
    if distt>max_deviation
         cell_register(i,2)=0;
    end
end

distt_mat=[];
for i=1:size(cell_register,1)
    if cell_register(i,2)>0
    distt_mat(i,1)=sum((c2(cell_register(i,2),:)-c1(i,:)).^2).^0.5;
    end
end