% Lopes-dos-Santos et al. 2013 modified ICA

function [icaSig,group]=Lopes_dos_Santos_2013(sig)

     sigz=zscore(sig,[],2);

     [coeff, score, latent, tsquared, explained, mu] = pca(sigz');

     % determine significant PCs
     [~,hbd,eigVals]=Marcenko_Pastur(score);
     numAssemb=sum(eigVals>hbd);

     pcaSig=score(:,1: numAssemb);

     [icaSig,A,W]=fastica(pcaSig');
     % find group by correlation
     group=[];
     for j=1:size(sigz,1)
          curSig=sigz(j,:);
          all_corr=zeros(size(icaSig,1),1);
          for i=1:size(icaSig,1)
             c=corrcoef(icaSig(i,:)',curSig');
             all_corr(i,1)=c(2);
          end
          if max(all_corr)>0.3 && max(all_corr)>quantile(all_corr(all_corr~=max(all_corr)),0.95)
              group(j)=find(all_corr==max(all_corr));
          else
              group(j)=-1;
          end
     end
 
 close all;