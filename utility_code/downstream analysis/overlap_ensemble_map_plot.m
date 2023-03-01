function rcomb=overlap_ensemble_map_plot(behav,group,firingrateAll_ensem_all,countTime_ensem_all)


    colorclust=distinguishable_colors(10);
    colorclust(4,:)=colorclust(10,:);

    rcomb=zeros(length(0:binsize:ceil(max(behav.position(:,2)))),length(0:binsize:ceil(max(behav.position(:,1)))),3);
    uni_gp=unique(group);


    for i=1:length(uni_gp)
    % 
        rs=filter2DMatrices(firingrateAll_ensem_all{i},1);
        r1=rs>0.5*max(rs(:));
        [ridxi,ridxj]=find(r1==1);

        for k1=1:length(ridxi)
            if sum(rcomb(ridxi(k1),ridxj(k1),:))==0
                rcomb(ridxi(k1),ridxj(k1),:)=colorclust(i,:);
            else
                rcomb(ridxi(k1),ridxj(k1),:)=(squeeze(rcomb(ridxi(k1),ridxj(k1),:))'+colorclust(i,:))./2;
            end

        end

        for k1=1:size(rcomb,3)
            rt=squeeze(rcomb(:,:,k1));
            rt(countTime_ensem_all==0)=255;
            rcomb(:,:,k1)=rt;
        end
    end
    
    imagesc(rcomb);
    axis ij
    shading flat;
    axis image
    axis off

