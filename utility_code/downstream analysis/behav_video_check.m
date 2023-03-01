function behav_video_check(behav,behav_prefix)

cd(behav.oripath);

fdir=dir([behav_prefix,'*.avi']);
for i=1:length(fdir)
    fname{i}=fdir(i).name;
end
fname=natsort(fname);
ctt=1;
v_all={};
tic;
for i=1:length(fdir)
    V=general_avi_read(fname{i});
    for j=1:length(V)
        v_all{ctt}=V{j};
        ctt=ctt+1;
    end
    toc;
end

pos_vid=behav.position;
pos_vid(:,1)=pos_vid(:,1)+behav.ROI(1);
pos_vid(:,2)=pos_vid(:,2)+behav.ROI(2);
pos_vid=pos_vid*behav.ROI3/behav.trackLength;

for i=1:length(v_all)
    imagesc(v_all{i});
    hold on;
    plot(pos_vid(1:i,1),pos_vid(1:i,2),'-','lineWidth',2);
    drawnow;
    clf
end

