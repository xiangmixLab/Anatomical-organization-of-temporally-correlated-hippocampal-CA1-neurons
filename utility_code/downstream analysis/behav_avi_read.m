function V=behav_avi_read(behavVidName)

V={};
ctt=1;
for i=1:length(behavVidName)
    Vtemp=general_avi_read(behavVidName{i});
    for j=1:length(Vtemp)
        V{ctt}=Vtemp{j};
        ctt=ctt+1;
    end
end
