% this cluster overlap is still based on rand index, but follows modi et
% al. 2014.
% original rand index = (a+b)/nchoosek(n,2), where a is the number of pairs
% stay together in both and b is the # of pairs separate in both
% This overlap is just (a)/nchoosek(n,2). no separate pairs considered
% 
function [overlap_idx_avgClust,overlap_idx,overlap_idx_perClust]=new_cluster_overlap_latest(c1,c2)

C=Contingency(c1,c2);	%form contingency matrix

%% total pairs of elements matched elements between c1 and c2
together_pairs=[];
for i=1:size(C,1)
    for j=1:size(C,2)
        if C(i,j)>2
            together_pairs(i,j)=nchoosek(C(i,j),2);
        end
        if C(i,j)<=2 && C(i,j)>0
            together_pairs(i,j)=1;
        end
        if C(i,j)==0
            together_pairs(i,j)=0; % if this item is 0, equals to Rand index's a
        end
    end
end

%% total pairs of c1, pre cluster
clust_pairs1=[];
for i=1:size(C,1)
    if sum(C(i,:))>2
        clust_pairs1(i,1)=nchoosek(sum(C(i,:)),2);
    end
    if sum(C(i,:))==2
        clust_pairs1(i,1)=1;
    end
    if sum(C(i,:))<2
        clust_pairs1(i,1)=1; % only one element in the cluster, set as 1 to avoid divide by 0
    end
end

%% total pairs of c2, pre cluster
clust_pairs2=[];
for i=1:size(C,2)
    if sum(C(:,i))>2
        clust_pairs2(i,1)=nchoosek(sum(C(:,i)),2);
    end
    if sum(C(:,i))==2
        clust_pairs2(i,1)=1;
    end
    if sum(C(:,i))<2
        clust_pairs2(i,1)=1; % only one element in the cluster, set as 1 to avoid divide by 0
    end
end

%% overlap calculation
overlap_idx_avgClust1=mean(sum(together_pairs,2)./clust_pairs1); % overlap for each cluster, then take average
overlap_idx_total1=sum(together_pairs(:))/sum(clust_pairs1); % overlap calculated by all overlap pairs / all possible intra-clsuter pairs
overlap_idx_perClust1=sum(together_pairs,2)./clust_pairs1;

overlap_idx_avgClust2=mean(sum(together_pairs,1)./clust_pairs2'); % overlap for each cluster, then take average
overlap_idx_total2=sum(together_pairs(:))/sum(clust_pairs2); % overlap calculated by all overlap pairs / all possible intra-clsuter pairs
overlap_idx_perClust2=sum(together_pairs,1)./clust_pairs2';

overlap_idx_avgClust=overlap_idx_avgClust1;
overlap_idx=overlap_idx_total1;
overlap_idx_perClust=overlap_idx_perClust1;

function Cont=Contingency(Mem1,Mem2)

if nargin < 2 | min(size(Mem1)) > 1 | min(size(Mem2)) > 1
   error('Contingency: Requires two vector arguments')
   return
end

Cont=zeros(max(Mem1),max(Mem2));

for i = 1:length(Mem1);
    if Mem1(i)>0&&Mem2(i)>0
        Cont(Mem1(i),Mem2(i))=Cont(Mem1(i),Mem2(i))+1;
    end
end
