% this cluster overlap is still based on rand index, but follows modi et
% al. 2014.
% original rand index = (a+b)/nchoosek(n,2), where a is the number of pairs
% stay together in both and b is the # of pairs separate in both
% This overlap is just (a)/nchoosek(n,2). no separate pairs considered
% 
function overlap_idx=new_cluster_overlap_032121(c1,c2)

C=Contingency(c1,c2);	%form contingency matrix

together_pairs=[];
for i=1:size(C,1)
    for j=1:size(C,2)
        if C(i,j)>2
            together_pairs(i,j)=nchoosek(C(i,j),2);
        end
        if C(i,j)>0&&C(i,j)<=2
            together_pairs(i,j)=1;
        end
        if C(i,j)==0
            together_pairs(i,j)=0;
        end
    end
end

for i=1:size(C,1)
    clust_pairs(i,1)=nchoosek(sum(C(i,:)),2); % total number of possible pairs in c1 cluster
end
idx_choose=sum(C,2)>1; % if there's a 1 element cluster, discard it from overlap calculation as the overlap will be 100%
sum_together_pairs=sum(together_pairs,2);
overlap_idx=mean(sum_together_pairs(idx_choose)./clust_pairs(idx_choose)); % overlap for each cluster, then take average (sum(together_pairs,2) indicates the number of possible pairs of c1 elements after distributed in c2



function Cont=Contingency(Mem1,Mem2)

if nargin < 2 | min(size(Mem1)) > 1 | min(size(Mem2)) > 1
   error('Contingency: Requires two vector arguments')
   return
end

Cont=zeros(max(Mem1),max(Mem2));

for i = 1:length(Mem1);
   Cont(Mem1(i),Mem2(i))=Cont(Mem1(i),Mem2(i))+1;
end
