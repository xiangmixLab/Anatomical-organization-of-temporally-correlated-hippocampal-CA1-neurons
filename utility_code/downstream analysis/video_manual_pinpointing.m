function [coor]=video_manual_pinpointing(vid,frameLabel)

coor=[]; % output will be x,y format
for i=frameLabel
    imagesc(vid(:,:,i));
    [x,y]=ginput
    x=mean(x);
    y=mean(y);
    coor=[coor;[x,y]];
end

coor=round(coor);