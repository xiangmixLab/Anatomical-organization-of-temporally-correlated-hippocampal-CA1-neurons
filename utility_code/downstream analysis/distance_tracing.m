function [newList,newListIdx]=distance_tracing(pointListOri)
    
    % part1: find seed
    seedIdx=1;
    for i=1:length(pointListOri)
        dist2item=sum((pointListOri-pointListOri(i,:)).^2,2).^0.5;
        if sum(dist2item<2)==2 % itself and the only neighbor
            seedIdx=i;
            break;
        end
    end

    
    %part 2:  for each point in pointList, reorder them in newList so that each point are next to each other in order (determined by seed selection)
    pointList=pointListOri;
    
    seed=pointList(seedIdx,:);

    pointList(seedIdx,:)=[];

    newList=[];
    newList(1,:)=seed;

    newListIdx=[];
    newListIdx(1,:)=seedIdx;

    ctt_newList=2;
    while ~isempty(pointList)

        dist2seed=sum((pointList-seed).^2,2).^0.5;
        
        min_value = min(dist2seed);
        
        if min_value<10
            min_index = find(dist2seed==min_value); % if multiple min_value exists, the first one is here
            min_index =min(min_index);
            newList(ctt_newList,:)=pointList(min_index,:);
            newListIdx(ctt_newList,1)=find(ismember(pointListOri, pointList(min_index,:),'rows'));

            ctt_newList=ctt_newList+1;
            
            seed=pointList(min_index,:);
            pointList(min_index,:)=[];
        else
            break;
        end
    end
        
    % last element chk
    if sum((newList(end,:)-newList(end-1,:)).^2)^0.5>2
        newList(size(newList,1),:)=[];
        newListIdx(size(newListIdx,1),:)=[];
    end

