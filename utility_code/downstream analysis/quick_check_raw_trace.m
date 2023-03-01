function raw_trace=quick_check_raw_trace(Y,A)

footprint=reshape(A,240,376)>0;

raw_trace=[];

for i=1:size(Y,3)
    ft=double(squeeze(Y(:,:,i))).*double(footprint);
    raw_trace(i)=sum(ft(:));
end

raw_trace(raw_trace==0)=nan;
raw_trace=fillmissing(raw_trace,'linear');

figure;
plot(raw_trace);