function [acolor,acolorreg,avgregion,avgregionshuf]=suppl17_tuncat_region_sz_util(group_ori,neuron)

acolor={};
acolorreg={};
avgregion={};
avgregionshuf={};
for k=2:10
    [acolor{k},acolorreg{k},avgregion{k},~,~,~,avgregionshuf{k}]=DBSCAN_region_quantify_022422(group_ori{k},{neuron},[]);      
end