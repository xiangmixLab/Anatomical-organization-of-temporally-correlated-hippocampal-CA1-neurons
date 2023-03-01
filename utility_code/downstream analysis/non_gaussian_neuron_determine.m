function indices=non_gaussian_neuron_determine(neuron,KL_thresh)

nA=neuron.A;
dshape=neuron.imageSize;
indices=[];

for i=1:size(neuron.A,2)
    [centroid,cov]=calculateCentroid_and_Covariance(nA(:,i),dshape(1),dshape(2));

    if ~isempty(centroid)&&~isempty(cov)&&(sum(sum(isnan(centroid)))==0||sum(sum(isnan(cov)))==0)
        A=mvnpdf_cal(dshape,centroid,cov);
        A1=nA(:,i);
        A=reshape(A,1,size(A,1)*size(A,2))';

        [KLt]=KLDiv_ori(A,A1);

        if KLt>KL_thresh
            indices(i)=1;
        end
    else
        indices(i)=1;
    end
end