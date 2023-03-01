
% adapted from Moser 2019 MEC objvec paper
% (https://www.nature.com/articles/s41586-019-1077-7), check the
% correlation of the deg-dis map in training and testing

% obj_tr and obj_ts corresponds to the same obj of different trials
function [objvec_cell,corr_val,shuf_corr_val]=objvec_cell_identification(fr_tr,fr_ts,obj_tr,obj_ts,binsize,corr_thresh)

    if ~isempty(fr_tr)&&~isempty(fr_ts)
        %% obj pos
        obj_tr=round(obj_tr/binsize);
        obj_ts=round(obj_ts/binsize);

        obj_tr(:,2)=size(fr_tr,1)-obj_tr(:,2);
        obj_ts(:,2)=size(fr_ts,1)-obj_ts(:,2);

        %% distance-orientation map
        [X_tr,Y_tr]=meshgrid(1:size(fr_tr,2),1:size(fr_tr,1));
        tr_pixlist=[reshape(X_tr,[],1),reshape(Y_tr,[],1)];

        diff2obj2tr=tr_pixlist-obj_tr;
        distance2obj_tr=round(sum((diff2obj2tr).^2,2).^0.5)+1;
        orientation2obj_tr=round(rad2deg(atan2(diff2obj2tr(:,1),diff2obj2tr(:,2)))+180)+1;

        [X_ts,Y_ts]=meshgrid(1:size(fr_ts,2),1:size(fr_ts,1));
        ts_pixlist=[reshape(X_ts,[],1),reshape(Y_ts,[],1)];

        diff2obj2ts=ts_pixlist-obj_ts;
        distance2obj_ts=round(sum((diff2obj2ts).^2,2).^0.5)+1;
        orientation2obj_ts=round(rad2deg(atan2(diff2obj2ts(:,1),diff2obj2ts(:,2)))+180)+1;

        dis_ori_map_tr=zeros(max(max(distance2obj_tr),max(distance2obj_ts)),361);
        dis_ori_map_ts=dis_ori_map_tr;

        fr_tr=filter2DMatrices(fr_tr,1);
        fr_ts=filter2DMatrices(fr_ts,1);

        fr_tr(fr_tr<0.5*max(fr_tr(:)))=0;
        fr_ts(fr_ts<0.5*max(fr_ts(:)))=0;
        
        for i=1:size(distance2obj_tr,1)
            dis_ori_map_tr(distance2obj_tr(i,1),orientation2obj_tr(i,1))=fr_tr(tr_pixlist(i,2),tr_pixlist(i,1));
        end
        for i=1:size(distance2obj_ts,1)
            dis_ori_map_ts(distance2obj_ts(i,1),orientation2obj_ts(i,1))=fr_ts(ts_pixlist(i,2),ts_pixlist(i,1));
        end

        %% corr
        dis_ori_map_tr1=filter2DMatrices(dis_ori_map_tr,1);
        dis_ori_map_ts1=filter2DMatrices(dis_ori_map_ts,1);

        dis_ori_map_tr1=reshape(dis_ori_map_tr1,[],1);
        dis_ori_map_ts1=reshape(dis_ori_map_ts1,[],1);

    %     dis_ori_map_tr1(dis_ori_map_tr1==0)=[];
    %     dis_ori_map_ts1(dis_ori_map_ts1==0)=[];

    %     dis_ori_map_ts1=resample(dis_ori_map_ts1,length(dis_ori_map_tr1),length(dis_ori_map_ts1));
    %     % we do not need to trim 0s, they represent an essential part of the
    %     map
        corrt=corrcoef(dis_ori_map_tr1,dis_ori_map_ts1);

        corr_val=corrt(2);

        %% shuffle vals
        shuf_corr_val=zeros(1,1000);
        for i=1:1000
            fr_tr_s=reshape(fr_tr(randperm(size(fr_tr,1)*size(fr_tr,2))),size(fr_tr,1),size(fr_tr,2));
            fr_ts_s=reshape(fr_ts(randperm(size(fr_ts,1)*size(fr_ts,2))),size(fr_ts,1),size(fr_ts,2));

            for j=1:size(distance2obj_tr,1)
                dis_ori_map_tr(distance2obj_tr(j,1),orientation2obj_tr(j,1))=fr_tr_s(tr_pixlist(j,2),tr_pixlist(j,1));
            end
            for j=1:size(distance2obj_ts,1)
                dis_ori_map_ts(distance2obj_ts(j,1),orientation2obj_ts(j,1))=fr_ts_s(ts_pixlist(j,2),ts_pixlist(j,1));
            end

            dis_ori_map_tr1=filter2DMatrices(dis_ori_map_tr,1);
            dis_ori_map_ts1=filter2DMatrices(dis_ori_map_ts,1);

            dis_ori_map_tr1=reshape(dis_ori_map_tr,[],1);
            dis_ori_map_ts1=reshape(dis_ori_map_ts,[],1);

    %         dis_ori_map_tr1(dis_ori_map_tr1==0)=[];
    %         dis_ori_map_ts1(dis_ori_map_ts1==0)=[];
    % 
    %         dis_ori_map_ts1=resample(dis_ori_map_ts1,length(dis_ori_map_tr1),length(dis_ori_map_ts1));
            corrt=corrcoef(dis_ori_map_tr1,dis_ori_map_ts1);

            shuf_corr_val(i)=corrt(2);
        end

        %% objvec determine
        if corr_val>quantile(shuf_corr_val,0.99)&&corr_val>corr_thresh
            objvec_cell=1;
        else
            objvec_cell=0;
        end
    else
        objvec_cell=0;
        shuf_corr_val=zeros(1,1000);
        corr_val=0;
    end