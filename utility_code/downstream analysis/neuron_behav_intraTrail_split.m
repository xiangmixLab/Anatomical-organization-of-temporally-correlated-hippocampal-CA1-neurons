
function [nnew,behavnew]=neuron_behav_intraTrail_split(neuron,behav,splitIdx)

% 1. original neuron and behav length
neuronLeng=size(neuron.C,2);
neuronTimeLeng=length(neuron.time);
behavLeng=size(behav.position,1);

% 2: neuron, behav resample to splitIdx
nC=[];
nCraw=[];

for i=1:size(neuron.C,1)
    nC(i,:)=resample(neuron.C(i,:),length(splitIdx),neuronLeng);
    nCraw(i,:)=resample(neuron.C_raw(i,:),length(splitIdx),neuronLeng);
end
ntime=resample(neuron.time,length(splitIdx),neuronTimeLeng);

bp=[resample(behav.position(:,1),length(splitIdx),behavLeng),resample(behav.position(:,2),length(splitIdx),behavLeng)];
bt=resample(behav.time,length(splitIdx),behavLeng);

% 3. neuron split
nnew={};
behavnew={};
usplitIdx=unique(splitIdx);
for i=1:length(usplitIdx)
    idxt=splitIdx==usplitIdx(i);
    nC1=nC(:,idxt);
    nCraw1=nCraw(:,idxt);
    ntime1=ntime(idxt);
    bp1=bp(idxt,:);
    bt1=bt(idxt,:);
    
    nnewTemp.C=[];
    nnewTemp.C_raw=[];
    nnewTemp.time=[];
    
    bnewTemp.position=[];
    bnewTemp.time=[];
    bnewTemp.ROI=behav.ROI;
    bnewTemp.trackLength=behav.trackLength;
    % 4: resample back to ori length
    for j=1:size(nC1,1)
        nnewTemp.C(j,:)=resample(nC1(j,:),neuronLeng,length(splitIdx));
        nnewTemp.C_raw(j,:)=resample(nCraw1(j,:),neuronLeng,length(splitIdx));
    end
    nnewTemp.S=C_to_peakS(nnewTemp.C);
    nnewTemp.time=resample(ntime1,neuronTimeLeng,length(splitIdx));

    bnewTemp.position=[resample(bp1(:,1),behavLeng,length(splitIdx)),resample(bp1(:,2),behavLeng,length(splitIdx))];
    bnewTemp.time=resample(bt1,behavLeng,length(splitIdx));
    
    nnewTemp.A=neuron.A;
    nnew{i}=nnewTemp;
    behavnew{i}=bnewTemp;
end

    

