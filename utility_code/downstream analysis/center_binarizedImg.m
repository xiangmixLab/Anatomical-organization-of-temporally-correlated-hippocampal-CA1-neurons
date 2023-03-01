function centeredImg=center_binarizedImg(binImg)

% 1. find the item
stats=regionprops(binImg);

% 2. find the offset
oriCen=[stats.BoundingBox(1)+stats.BoundingBox(3)/2,stats.BoundingBox(2)+stats.BoundingBox(4)/2];
realCen=size(binImg)/2;
offset=[realCen(2)-oriCen(1),realCen(1)-oriCen(2)];

% 3. generate centeredImg
centeredImg=zeros(size(binImg));
[i,j]=find(binImg==1);
centeredImg(i+offset(2),j+offset(1))=1;

centeredImg=logical(centeredImg);