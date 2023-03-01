function cidx_str=num2str_largeNum(cidx)
cidx_str=[];
while 1
    cidx1=mod(cidx,10);
    cidx_str=[num2str(cidx1),cidx_str];
    cidx=floor(cidx/10);
    if cidx<1
        break;
    end
end