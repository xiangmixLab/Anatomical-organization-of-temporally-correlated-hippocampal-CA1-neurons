%% Dombeck & Tank 2010 2-p VR place cell identification paper method (they use dF/F, substitute to firing rate)

function [place_cells,firingrate]=place_cell_identification_tank2010(neuron,behavpos,behavtime,ROI,binsize,temp,small_velo)

    [firingrate,count,~,countTime] = calculatingCellSpatialForSingleData_Suoqin(neuron,behavpos,behavtime,ROI,binsize,1:size(neuron.C,1),0.1*max(neuron.C,[],2),temp,[],[],[0.1 inf],small_velo);
    
    overall_mean_rate=zeros(length(count),1);
    for i=1:length(count)
        overall_mean_rate(i)=nansum(count{i}(:))/nansum(countTime(:));
    end
    
    time0=neuron.time;
    
    rng default % for reproduction

    firingrate=filter2DMatrices_021521(firingrate,5,2);   
    [~,numFields]=place_field_define(firingrate);
    
    % 6: survive the shuffling
    % rotate_leng=floor(size(S0,2)/(nboot+1));
    C_dat=trunk_split_data(neuron.C,100,'col');
    S_dat=trunk_split_data(C_to_peakS(neuron.C),100,'col');
    
    for r=1:100
        neuronboot=trunk_shuffle_data_neuron_simple(C_dat,S_dat,time0);
        [firingrateAllt{r}] = calculatingCellSpatialForSingleData_Suoqin(neuronboot,behavpos,behavtime,ROI,binsize,1:size(neuronboot.C,1),0.1*max(neuronboot.C,[],2),temp,[],[],[0.1 inf],small_velo);       
        firingrateAllt{r}=filter2DMatrices_021521(firingrateAllt{r},5,2);
    end
    
    numFieldsBoot=[];
    for r=1:100         
        [~,numFieldsBoot(:,r)]=place_field_define(firingrateAllt{r});
    end
        
    place_cells=quantile(numFieldsBoot,0.95,2)==0;        
    place_cells(numFields==0)=0;
    place_cells(logical(double(numFields==0).*double(overall_mean_rate<0.01)))=0;
    place_cells=find(place_cells==1);
end

function [pf_all,numFields_all,numFields_ori]=place_field_define(fr_all)
    
    pf_all=cell(length(fr_all),1);
    numFields_all=zeros(length(fr_all),1);
    
    for i=1:length(fr_all)
        
        fr=fr_all{i};
        
        if ~isempty(fr)&&sum(fr(:))>0
            % 1: bin rate should be larger than 0.25*(max(fr)-percentile(fr,lower 0.25))
            fr_b=fr>max(fr(:))*0.25;

            % 2: place field should larger than 3*3 bins (2-3cm*2-3cm) and smaller
            % than 1/9 of total size
            potentialFields={};

            statss=regionprops(fr_b,'PixelIdxList','BoundingBox','Area');
            ctt=1;
            for j=1:length(statss)
                if statss(j).BoundingBox(3)>=3&&statss(j).BoundingBox(4)>=3 && statss(j).BoundingBox(3)<size(fr_b,2)/2&&statss(j).BoundingBox(4)<size(fr_b,1)/2
                    potentialFields{ctt}=statss(j);
                    ctt=ctt+1;
                end
            end

            del_idx=zeros(length(potentialFields),1);

            % 3: mean field rate should be larger than 3*out_field_rate
            for j=1:length(potentialFields)
                fr_field_mean=mean(mean(fr(potentialFields{j}.PixelIdxList)));
                if fr_field_mean<=mean(mean(fr(setdiff([1:size(fr,1)*size(fr,2)],potentialFields{j}.PixelIdxList))))*2
                    del_idx(j)=1;
                end
            end

            numFields_ori=length(potentialFields);
            potentialFields(logical(del_idx))=[];
            numFields=length(potentialFields);
            pf_all{i}=potentialFields;
            numFields_all(i)=numFields;
            
        else
            pf_all{i}=[];
            numFields_all(i)=0;          
        end
    end
end