function [groupAnatomical,patches]=major_anatomical_clusters(group1,centroidd)

gp1=group1(group1~=-1);
anatomicalGpCounter=11;

groupAnatomical=group1*0;
groupAnatomical(group1==-1)=-1;

neighboor_dis_all=squareform(pdist(centroidd));
neighboor_dis_all(neighboor_dis_all==0)=inf;
neighboor_dis=quantile(min(neighboor_dis_all,[],2),0.95);

for gp=1:length(unique(gp1))

    cen=round(centroidd(group1==gp,:));
    nd=neighboor_dis;

    [IDX, isnoise]=DBSCAN_modified(cen,nd,4,'dis');
    
    uIDX=unique(IDX);
    uIDX(uIDX==0)=[];
    
    groupAnatomical(group1==gp)=IDX;
    for i=1:length(uIDX)
        groupAnatomical(groupAnatomical==uIDX(i))=anatomicalGpCounter;
        anatomicalGpCounter=anatomicalGpCounter+1;
    end
end

groupAnatomical=groupAnatomical-10;
groupAnatomical(groupAnatomical<=0)=-1;

% major anatomical cluster patches
patches={};
for gp=1:length(unique(gp1))
    cen=round(centroidd(group1==gp,:));
    anatomical_region_idx=groupAnatomical(group1==gp);
    u_region_idx=unique(anatomical_region_idx);
    for j=1:length(u_region_idx)
        region_cen=cen(anatomical_region_idx==u_region_idx(j),:);
        region_bd=region_cen(boundary(region_cen(:,1),region_cen(:,2)),:);
        patches{gp,j}=region_bd;
    end
end