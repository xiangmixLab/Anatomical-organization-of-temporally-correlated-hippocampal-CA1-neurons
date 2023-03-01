function [possibilities,total_combination]=total_possible_cluster_organizations(group)

uni_group=unique(group);
total_elements=length(group);

num_group=[];
for i=1:length(uni_group)
    num_group(i)=sum(group==uni_group(i));
end

total_combination=1;
for i=1:length(uni_group)
    total_combination=total_combination*nchoosek(total_elements-sum(num_group(1:i-1)),num_group(i))/10^69; 
    % 1: first group: suppose it has K element, total num of elements is N, so the total possible pickups of first group is C1=(C^K_N).
    % 2: next group: suppose it has K1 element, total num of elements is N-K,so the total possible pickups of first group is C2=(C^K1_(N-K)).
    % the total possible pickups of the two group: C1*C2
end

% validation: change calculation direction
total_combination1=1;
for i=length(uni_group):-1:1
    total_combination1=total_combination1*nchoosek(total_elements-sum(num_group(i+1:length(uni_group))),num_group(i))/10^50; 
    % 1: first group: suppose it has K element, total num of elements is N, so the total possible pickups of first group is C1=(C^K_N).
    % 2: next group: suppose it has K1 element, total num of elements is N-K,so the total possible pickups of first group is C2=(C^K1_(N-K)).
    % the total possible pickups of the two group: C1*C2
end    

% possibilities
possibilities=1/total_combination; % P(clust)