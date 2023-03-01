function  [PF_radius,other_fields,other_field_stats,Areaa,localmax,autocorr]=firing_field_radius_autocorr(fr,smooth_sign)

do_sign=~isempty(fr);
if ~isempty(fr)&&isnan(fr(1))
    do_sign=0;
end

if do_sign
    if smooth_sign==1
        fr1=filter2DMatrices(fr,1); % fr will be non-smoothed upon request
    else
        fr1=fr;
    end
    autocorr=xcorr2(fr1);
    % PF_radius= findPlaceFieldRadius(autocorr);
    autocorr_bin=autocorr>=0.5*max(autocorr(:));
    
    stats=regionprops(autocorr_bin,'Area','MajorAxisLength');

    all_area=[stats.Area];
    if min(size(fr))>1
        all_dia=[stats.MajorAxisLength]; % half width at half height of original autocorr plot
    else
        all_dia=all_area;
    end
        
    [all_area_sort,idx]=sort(all_area); % if have multiple fields
    idx_large_fields=[1, find(all_area_sort(2:end)>all_area_sort(1)*0.8)+1]; % find those with competetive field size with the largest one

%     all_dia_sort=all_dia(idx);
    PF_radius=mean(all_dia)/2;

    all_area_select=all_area_sort(idx_large_fields);
    Areaa=mean(all_area_select);

    autocorr_remain=autocorr;
    autocorr_remain(autocorr>=0.5*max(autocorr(:)))=0.5*max(autocorr(:));
    other_fields=imregionalmax(autocorr_remain);
    other_field_stats=regionprops(other_fields);

    high_field_stats=regionprops(autocorr_bin);
    for k=1:length(other_field_stats)
        for j=1:length(high_field_stats)
            if norm([other_field_stats(k).Centroid;high_field_stats(j).Centroid])<1
                break;
            end
        end
    end

    other_field_stats(k)=[];
    
    localmaxMap=imregionalmax(autocorr);
    stats1=regionprops(localmaxMap);
    localmax=length(stats1);
else
    autocorr=[];
    PF_radius=nan;
    other_fields=[];
    other_field_stats=[];
    localmax=nan;
    Areaa=nan;
end