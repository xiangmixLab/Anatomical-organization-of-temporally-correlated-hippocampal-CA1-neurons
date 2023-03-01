function [KL]=neuron_KL_calculation(neuron)

nA=neuron.A;
dshape=neuron.imageSize;
KL=[];

parfor i=1:size(neuron.A,2)
    [centroid,cov]=calculateCentroid_and_Covariance(nA(:,i),dshape(1),dshape(2));

    if ~isempty(centroid)&&~isempty(cov)&&(sum(sum(isnan(centroid)))==0||sum(sum(isnan(cov)))==0)
        A=mvnpdf_cal(dshape,centroid,cov);
        A1=nA(:,i);
        A=reshape(A,1,size(A,1)*size(A,2))';

        [KLt]=KLDiv_ori(A,A1);
        KL(i)=KLt;
    end
end