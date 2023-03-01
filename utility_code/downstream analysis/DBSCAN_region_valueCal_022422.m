function [avg_region,max_region,avg_region_nnum,max_region_nnum,max_regions,area_thresh_all]=DBSCAN_region_valueCal_022422(region_all,region_nnum_all,imageSize)

d1=imageSize(1);
d2=imageSize(2);

avg_area=[];
max_area=[];
avg_area_nnum=[];
max_area_nnum=[];
area_thresh_t=[];

for k=1:length(region_all)
    regiont_all=zeros(d1,d2);
    for m=1:length(region_all{k})
        if ~isempty(region_all{k}{m})
            regiont_all=double(regiont_all)+double(region_all{k}{m});
        end
    end
    regiont_all=regiont_all>0;
    stats=regionprops(regiont_all,'Area','Perimeter');
    
    region_nnum_t=cell2mat(region_nnum_all{k});
    
    if ~isempty(stats)
        area_all=[stats.Area];
        area_all_thres=max(area_all)*0.1; % 95% max 122019 % change back to 0.1 012720
        area_idx=area_all>area_all_thres;
        area_thresh_t(k)=area_all_thres;
        avg_area(k)=mean(area_all(area_idx));
        max_area(k)=max(area_all(area_idx));

        num_thres=max(region_nnum_t)*0.1;
        avg_area_nnum(k)=ceil(mean(region_nnum_t(region_nnum_t>num_thres)));
        max_area_nnum(k)=ceil(max(region_nnum_t(region_nnum_t>num_thres)));
    else
        avg_area(k)=1;
        avg_area_nnum(k)=0;
        max_area_nnum(k)=0;
        max_area(k)=1;
        area_thresh_t(k)=0; % no regions actually...
    end
end
area_thresh_t(isnan(area_thresh_t))=0;
area_thresh_t(area_thresh_t==inf)=0;
area_thresh_all=area_thresh_t;

avg_region=mean(avg_area);
max_region=mean(max_area);
avg_region_nnum=mean(avg_area_nnum);
max_region_nnum=mean(max_area_nnum);
max_regions={avg_area,max_area};

