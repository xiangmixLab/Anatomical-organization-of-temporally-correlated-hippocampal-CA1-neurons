function [velo,ddistance,dtime,cen,diff_from_cen]=behav_velo_cal(behav,cen,led)

behav.position(:,1)=fillmissing(behav.position(:,1),'nearest');
behav.position(:,2)=fillmissing(behav.position(:,2),'nearest');
behav.time=fillmissing(behav.time,'linear');

if isequal(led,'r')
    for i=2:size(behav.position,1)
        ddistance(i-1,1)=norm(behav.position(i-1,:)-behav.position(i,:));
    end    
    dtime=diff(behav.time);    
    velo=ddistance./(dtime/1000); % convert to s
    if isempty(cen)
        cen=[mean(behav.position(:,1)),mean(behav.position(:,2))];
        for i=1:size(behav.position,1)
            diff_from_cen(i)=norm(behav.position(i,:)-cen);
        end
    else
        for i=1:size(behav.position,1)
            diff_from_cen(i)=norm(behav.position(i,:)-cen);
        end
    end
end

if isequal(led,'b')
    for i=2:size(behav.positionblue,1)
        ddistance(i-1,1)=norm(behav.positionblue(i-1,:)-behav.positionblue(i,:));
    end    
    dtime=diff(behav.time);    
    velo=ddistance./(dtime/1000);
    if isempty(cen)
        cen=[mean(behav.positionblue(:,1)),mean(behav.positionblue(:,2))];
        for i=1:size(behav.positionblue,1)
            diff_from_cen(i)=norm(behav.positionblue(i,:)-cen);
        end
    else
        for i=1:size(behav.positionblue,1)
            diff_from_cen(i)=norm(behav.positionblue(i,:)-cen);
        end
    end
end

