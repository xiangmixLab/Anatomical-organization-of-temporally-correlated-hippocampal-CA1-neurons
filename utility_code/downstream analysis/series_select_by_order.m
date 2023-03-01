% return ser1 elements appears in ser2 by order in ser1 (ser2 should be
% included in ser1)
function ser1_select=series_select_by_order(ser1,ser2)
ser1_select=[];
for i=1:length(ser1) 
    idxx=find(ser2==ser1(i));
    ser1_select=[ser1_select,ser2(idxx)];
end
        