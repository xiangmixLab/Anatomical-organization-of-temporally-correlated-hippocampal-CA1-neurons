function [avg_max_region,avg_max_peri,avg_max_region_95,avg_max_peri_95,area_thresh_all,all_mean_regions,all_max_regions,avg_max_region_nnum,all_nnums]=DBSCAN_region_perimeter_quantify_func_052821(region_all,region_nnum_all)

avg_max_region=[];
avg_max_region_nnum=[];
avg_max_peri=[];
area_thresh_all={};
all_mean_regions={};
all_max_regions={};

for i=1:size(region_all,1)
    for j=1:size(region_all,2)
        regiont=region_all{i,j};
        region_nnum_t=region_nnum_all{i,j};

        avg_area=[];
        avg_area_nnum=[];
        avg_peri=[];
        area_thresh_t=[];
        for k=1:size(regiont,1)
            area_all=[];
            area_all_nnum=[];
            peri_all=[];
             
            for l=1:size(regiont,2)
                if ~isempty(regiont{k,l})
                    regiont_all=regiont{k,l};
                    regiont_all=regiont_all>0;
                    stats=regionprops(regiont_all,'Area','Perimeter');
                    if ~isempty(stats)
                        area_all=[area_all,sum([stats.Area])];
                        area_all_nnum=[area_all_nnum,region_nnum_t{k,l}];
                        peri_all=[peri_all,sum([stats.Perimeter])];
                    end
                end
            end
            if ~isempty(area_all)
                area_all_thres=max(area_all)*0.2; % 95% max 122019 % change back to 0.1 012720 % change to 0.2 053021
                area_idx=area_all>area_all_thres;
                area_thresh_t(k)=area_all_thres;
                avg_area(k)=nanmean(area_all(area_idx));
                avg_area_nnum(k)=ceil(mean(area_all_nnum(area_idx)));
                max_area(k)=max(area_all(area_idx));
                max_area_nnum(k)=ceil(max(area_all_nnum(area_idx)));

                avg_peri(k)=mean(peri_all(area_idx));    
            else
                avg_area(k)=nan;
                avg_area_nnum(k)=nan;
                avg_peri(k)=nan;
            end
        end

        area_thresh_t(isnan(area_thresh_t))=0;
        area_thresh_t(area_thresh_t==inf)=0;
        area_thresh_all{i,j}=area_thresh_t;
        avg_max_region(i,j)=nanmean(avg_area);
        avg_max_region_nnum(i,j)=nanmean(max_area_nnum);
        all_mean_regions{i,j}=avg_area;
        all_max_regions{i,j}=max_area;
        all_nnums{i,j}=avg_area_nnum;
        avg_max_region_95(i,j)=quantile(avg_area,.95);
        avg_max_peri(i,j)=nanmean(avg_peri);
        avg_max_peri_95(i,j)= quantile(avg_peri,.95);
    end
end
