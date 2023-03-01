function [sumFiringRate_conv,sumFiringRateObject,sumFiringRate] = comparingFiringRateSingleConditionMultiObjects_110420(firingRate,countTime,objpos,binsize,behavROI)
    
    posObjects = round(objpos./binsize);

    [sumFiringRate,sumFiringRate_conv,cellnum] = calculatingSumFiringRate(firingRate,countTime);

    sumFiringRate=sumFiringRate/cellnum;
    sumFiringRate_conv=sumFiringRate_conv/cellnum;
%     sumFiringRate_conv(sumFiringRate_conv==0)=nan;
    
    sumFiringRateObject = zeros(1,size(posObjects,1));

    if isempty(behavROI)
        bROI=[0 0 size(sumFiringRate,2) size(sumFiringRate,1)];
    else
        bROI=round(behavROI/binsize);
    end
    
    obj_range=1; % an estimate of the range red light position may indicate that the mouse is interacting with obj
    
    if sum(posObjects(:))>0
        for i = 1:size(posObjects,1)
    %         sumFiringRateObject(1,i) =sumFiringRateObject(1,i)+sumFiringRate_conv(size(sumFiringRate_conv,1)-posObjects(i,2)+1+u,posObjects(i,1)+v);
            sumFiringRateObject(1,i) =0;
            sum_num=0;
            for u=-obj_range:obj_range % the obj is large, a larger area may be required (5x5 around obj)
                for v=-obj_range:obj_range
                    if bROI(4)-posObjects(i,2)+u>0&&bROI(4)-posObjects(i,2)+u<size(sumFiringRate_conv,1)&&posObjects(i,1)+v>0&&posObjects(i,1)+v<size(sumFiringRate_conv,2)
                        if ~isnan(sumFiringRate_conv(bROI(4)-posObjects(i,2)+u,posObjects(i,1)+v))
                            sumFiringRateObject(1,i) =sumFiringRateObject(1,i)+sumFiringRate_conv(bROI(4)-posObjects(i,2)+u,posObjects(i,1)+v); % change from sunFiringRate to sumFiringRate_conv
                            sum_num=sum_num+1;
                        end
                    end
                end
            end
            sumFiringRateObject(1,i)=sumFiringRateObject(1,i)/sum_num;
        end
    end
%     sumFiringRate_conv(size(sumFiringRate_conv,1)-posObjects(1,2)+1,posObjects(1,1))=-1;
%     sumFiringRate_conv(size(sumFiringRate_conv,1)-posObjects(2,2)+1,posObjects(2,1))=-1;
    sumFiringRateObject(isnan(sumFiringRateObject)) = 0;

end


function [sumFiringRate,sumFiringRate_conv,cellnum] = calculatingSumFiringRate(firingRate,countTime)

    
    radius=3;
%     f=ones(2*radius-1); % this is sum
    f=ones(2*radius-1)/sum(sum(ones(2*radius-1))); % this is average
    
    sumFiringRate = zeros(size(countTime));
    sumFiringRate_conv=zeros(size(countTime));
    
    cellnum=0;

    for i = 1:length(firingRate)
            if isempty(firingRate{i})
                continue;
            end
            firingRate{i}(isnan(firingRate{i})) = nan;
            firingRate{i}(firingRate{i}==inf) = nan;
            firingRate{i}(countTime==0) = nan;
            firingRate{i}(isnan(firingRate{i}))=0;
            sumFiringRate_conv = sumFiringRate_conv + nanconv(firingRate{i},f);
            sumFiringRate = sumFiringRate + firingRate{i};
            cellnum=cellnum+1;
    end
end
    % sumFiringRate(countTime == 0) = NaN;



