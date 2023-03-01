function outdat=idx_dupliucate_interleaver(dat,rep)

outdat=zeros(length(dat)*rep,1);
for i=1:rep
    outdat(i:rep:length(dat)*rep-rep+i)=dat;
end