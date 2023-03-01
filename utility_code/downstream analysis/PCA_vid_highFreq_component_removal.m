function Y_denoised=PCA_vid_highFreq_component_removal(Y)

% Yt=imresize(Y,0.5);
% Yt=Yt(:,:,1:2:end);
Y1=reshape(Y,size(Y,1)*size(Y,2),size(Y,3));

[coeff,score,latent,tsquared,explained,mu] = pca(Y1');

coeff(:,4:end)=coeff(:,4:end)*0;

Y_denoised_t=score*coeff' + repmat(mu,size(Y1,2),1);

Y_denoised=reshape(Y_denoised_t,size(Y,1),size(Y,2),size(Y,3));