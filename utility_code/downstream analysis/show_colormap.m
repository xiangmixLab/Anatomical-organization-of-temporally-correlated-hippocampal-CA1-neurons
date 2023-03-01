% input: colorMap, n*3 matrix with each row represent a RGB color

function show_colormap(colorMap)

colorMapImg=zeros(size(colorMap,1),1,3);
for i=1:size(colorMap,1)
    colorMapImg(i,1,:)=colorMap(i,:);
end

imagesc(colorMapImg);