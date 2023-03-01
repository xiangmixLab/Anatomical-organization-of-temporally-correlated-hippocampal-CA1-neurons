function gp1=clusterIdx_swap(group,i,j)

gp1=group;
gp1(gp1==i)=100;
gp1(gp1==j)=i;
gp1(gp1==100)=j;