function cellDat1=empty_cell_clean(cellDat)

empIdx=[];
for i=1:length(cellDat)
    if isempty(cellDat)
        empIdx=[empIdx,i];
    end
end

cellDat(empIdx)=[];
cellDat1=cellDat;
