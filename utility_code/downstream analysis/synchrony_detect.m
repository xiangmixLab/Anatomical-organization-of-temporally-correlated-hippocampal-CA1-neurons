function all_synchrony=synchrony_detect(neuron,group,Fs)

all_synchrony={};
ugroup=unique(group);
ugroup(ugroup==-1)=[];

for i=1:length(ugroup)
    nC=neuron.C(group==ugroup(i),:);
    nS=neuron.S(group==ugroup(i),:);
    nC_peak=C_to_peakS(nC);
    
    thresh=3*std(nS,[],2);
    
    nC_peak(nC_peak<thresh)=0;
    
    nEvent=nC_peak;
    win=Fs*50; % 50sec window
    
    all_synchrony{i}=[];
    for j=1:size(nEvent,2)-win
        curSegment=nEvent(:,j:j+win-1);
        curFiringSequence=[];
        for j1=1:size(curSegment,1)
            if sum(curSegment(j1,:))>0
                curFiringSequence(j1,1)=median(find(curSegment(j1,:)==max(curSegment(j1,:))));
            else
                curFiringSequence(j1,1)=-1;
            end
        end
        c1=curFiringSequence;
        c1(c1==-1)=[];
        if max(c1)-min(c1)<=10*Fs % all max response in 10 sec, a seq identified
            all_synchrony{i}=[all_synchrony{i},curFiringSequence];
        end
    end
    
    
end