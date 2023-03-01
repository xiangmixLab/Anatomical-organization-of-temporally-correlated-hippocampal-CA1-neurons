% https://www.mathworks.com/matlabcentral/answers/223299-error-the-variable-in-a-parfor-cannot-be-classified
function [A_color,A_color_region,avg_region_LT_tun,avg_region_shuf_LT_tun]=experiment_cluster_region_cal_2_to_10(group_ori,dpath)

miceNum=length(dpath);
trialNum=size(group_ori,2);

avg_region_LT_tun=cell(miceNum,trialNum);
avg_region_shuf_LT_tun=cell(miceNum,trialNum);
A_color=cell(miceNum,trialNum);
A_color_region=cell(miceNum,trialNum);
all_neurons=cell(miceNum,trialNum);

t1=tic;
for j=1:miceNum
    load([dpath{j},'\','neuronIndividuals_new_tuncat_cleaned.mat'])
    for j1=1:trialNum
        all_neurons{j,j1}=Sources2D;
        all_neurons{j,j1}.C=neuronIndividuals_new_tuncat{j1}.C;
        all_neurons{j,j1}.S=neuronIndividuals_new_tuncat{j1}.S;
        all_neurons{j,j1}.A=neuronIndividuals_new_tuncat{j1}.A;
        all_neurons{j,j1}.imageSize=[240,376];
        avg_region_LT_tun{j,j1}={};
        avg_region_shuf_LT_tun{j,j1}={};
        A_color{j,j1}={};
        A_color_region{j,j1}={};
    end
end
disp(['finish loading neuron, ',num2str(toc(t1)),'sec'])

%% generate idx
idxs={};
ctt=1;
for j=1:miceNum
    for j1=1:trialNum     
        idxs{ctt}=[j,j1];
        ctt=ctt+1;
    end
end

%% parfor
nums=length(idxs);
t1=tic;
a={};
b={};
c={};
d={};
parfor j=1:nums
     idx_curr=idxs{j};
     [acolor,acolorreg,avgregion,~,~,~,avgregionshuf]=DBSCAN_region_quantify_022422(group_ori{idx_curr(1),idx_curr(2)},{all_neurons{idx_curr(1),idx_curr(2)}},[]);      

     a{j}=acolor;
     b{j}=acolorreg;
     c{j}=avgregion;
     d{j}=avgregionshuf;
end   
disp(['finish, ',num2str(toc(t1))]);

%% insert res
ctt=1;
for j=1:miceNum
    for j1=1:trialNum
        idx_curr=idxs{ctt};
        A_color{idx_curr(1),idx_curr(2)}=a{ctt};
        A_color_region{idx_curr(1),idx_curr(2)}=b{ctt};
        avg_region_LT_tun{idx_curr(1),idx_curr(2)}=c{ctt};
        avg_region_shuf_LT_tun{idx_curr(1),idx_curr(2)}=d{ctt};
        ctt=ctt+1;
    end
end