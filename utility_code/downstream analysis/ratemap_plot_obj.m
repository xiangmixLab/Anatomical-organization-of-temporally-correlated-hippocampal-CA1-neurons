function ratemap_plot_obj(fr,countTime,smoothsign,colorbar_sign,obj,binsize,behavROI)
if smoothsign==1
    fr=filter2DMatrices_021521(fr,5,2);
end
% fr=imresize(fr,size(countTime));
fr(countTime==0)=nan;

if isempty(behavROI)
    bROI=[0 0 size(fr,2) size(fr,1)];
else
    bROI=round(behavROI/binsize);
end

pcolor(fr);
hold on;
if ~isempty(obj)
    for i=1:size(obj,1)
        obj(i,:)=round(obj(i,:)./binsize);
        obj(i,2)=bROI(4)-obj(i,2); % axes problems
        plot(obj(i,1),obj(i,2),'.','MarkerSize',24,'color','k');
        hold on;
    end
end

if colorbar_sign
    colorbar;
end

colormap(jet);
axis ij
shading flat;
axis image
axis off

