function corrP=correlation_significance(nC)
ctt=1;
for i=1:size(nC,1)-1
    for j=i+1:size(nC,1)
        [~,p]=corrcoef(nC(i,:),nC(j,:));
        corrP(ctt)=p(2);
        ctt=ctt+1;
    end
end

corrP=squareform(corrP);