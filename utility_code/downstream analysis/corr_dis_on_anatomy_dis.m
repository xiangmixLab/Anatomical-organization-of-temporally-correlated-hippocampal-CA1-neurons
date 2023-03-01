function [close_all,far_all,neighboor_dis]=corr_dis_on_anatomy_dis(neuronIndividuals,session,anatomical_close_dis)

% distance

distt=(squareform(pdist(neuronIndividuals{session}.centroid)));
distt(distt==0)=nan;
neighboor_dis_all=nanmin(distt,[],1);
neighboor_dis=quantile(neighboor_dis_all,0.7);

if ~isempty(anatomical_close_dis)
    neighboor_dis=anatomical_close_dis;
end

% correlation
dataC1=neuronIndividuals{session}.C;
dataC1(isnan(dataC1(:)))=0;
dataC1=zscore(dataC1,[],2);

corr_all=squareform(1-pdist(dataC1,'correlation'));

close_all=corr_all(distt<=neighboor_dis);
far_all=corr_all(distt>neighboor_dis);
