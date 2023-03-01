function [group12,group21]=determineReorganizedCluster_072922(group1,group2)

negativeOne_list=group1*0;
negativeOne_list(group1==-1)=1;
negativeOne_list(group2==-1)=1;

group1(logical(negativeOne_list))=[];
group2(logical(negativeOne_list))=[];

uni_g1=unique(group1);
uni_g2=unique(group2);

% represent how a cluster in group1 is consituted in group2
% group12 column1: group1 index
% group12 column2: group2 indexs
% group12 column3: cell index for each of the group2 indexs in column2
group12={};
for i=1:length(uni_g1)
    idxx=find(group1==uni_g1(i));
    g2_idxx=group2(idxx);
    uni_g2_idxx=unique(g2_idxx);
    group12{i,1}=i;
    group12{i,2}=uni_g2_idxx;
    
    gtemp={};
    for j=1:length(uni_g2_idxx)
        gtemp{j}=idxx(g2_idxx==uni_g2_idxx(j));
    end
    group12{i,3}=gtemp;
end

% represent how a cluster in group2 is consituted in group1
% group21 column1: group2 index
% group21 column2: group1 indexs
% group21 column3: cell index for each of the group1 indexs in column2

group21={};
for i=1:length(uni_g2)
    idxx=find(group2==uni_g2(i));
    g1_idxx=group1(idxx);
    uni_g1_idxx=unique(g1_idxx);
    group21{i,1}=i;
    group21{i,2}=uni_g1_idxx;
    
    gtemp={};
    for j=1:length(uni_g1_idxx)
        gtemp{j}=idxx(g1_idxx==uni_g1_idxx(j));
    end
    group21{i,3}=gtemp;
end