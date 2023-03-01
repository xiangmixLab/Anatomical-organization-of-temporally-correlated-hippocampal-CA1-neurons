function nS=C_to_peakS_with_thresh(nC,thresh)

nS=[];
for i=1:size(nC,1)
    t=nC(i,:);
    [pks,loc]=findpeaks(t);
    t=t*0;
    t(loc)=pks;
    t(t<thresh)=0;
    nS(i,:)=t;
end