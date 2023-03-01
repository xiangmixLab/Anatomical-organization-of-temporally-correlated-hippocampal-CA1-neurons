function position=behav_locate_blackmice(V,behav)

% suppose V is the full vid

% 1. calculate bg
bg=double(V{1});
for i=2:length(V)
    bg=bg+double(V{i});
end

bg=(bg/length(V));

% 2. substract bg
V1=V;
for i=1:length(V)
    V1{i}=double(V1{i})-(bg);
    V1{i}(V1{i}>0)=0;
    V1{i}=0-V1{i};
    V1{i}=imadjustRGB2017(V1{i});
end

% 3. find highlighted mice
ROI=behav.ROI/behav.trackLength*behav.ROI3;
se=strel('disk',17);
miceMask={};
for i=1:length(V1)
    rgbframe=double(V1{i});

    miceMaskt=double(rgbframe(:,:,1)>150).*double(rgbframe(:,:,2)>150).*double(rgbframe(:,:,3)>150);

    miceMask{i}=bwareaopen(miceMaskt,20);
    miceMask{i}=imclose(miceMask{i},se);
    miceMask{i}=bwareaopen(miceMask{i},130);
    miceMaskROI=miceMask{i}*0;
    miceMaskROI(ROI(2)+16:ROI(2)+ROI(4)-16,ROI(1)+16:ROI(1)+ROI(3)-16)=1;
    miceMask{i}=miceMask{i}.*miceMaskROI;
end
miceMask{1}=miceMask{1}*0;

% 4. Karman tracking
miceMaskV=concatenate_img_cell(miceMask,0);
alltracks_mouse = MotionBasedMultiObjectTrackingExample_adapted(miceMaskV,[]);

%% 5. arrange tracking ids
% total length
tic;
ids=[];
for i=1:length(alltracks_mouse)
    for j=1:length(alltracks_mouse{i})
        ids=[unique(ids);alltracks_mouse{i}(j).id];
    end
end
toc;
trackingRes=cell(max(ids),1);
trackingFrameStamp=cell(max(ids),1);

% assign track
for i=1:length(alltracks_mouse)
    for j=1:length(alltracks_mouse{i})
        coor=[alltracks_mouse{i}(j).bbox(1)+alltracks_mouse{i}(j).bbox(3)/2,alltracks_mouse{i}(j).bbox(2)+alltracks_mouse{i}(j).bbox(4)/2];
        if coor(1)>ROI(1)&&coor(1)<ROI(1)+ROI(3)&&coor(2)>ROI(2)&&coor(2)<ROI(2)+ROI(4)
            trackingRes{alltracks_mouse{i}(j).id}=[trackingRes{alltracks_mouse{i}(j).id};coor];
            trackingFrameStamp{alltracks_mouse{i}(j).id}=[trackingFrameStamp{alltracks_mouse{i}(j).id};i];
        end
    end
end

%% 6. find robust tracks
position=nan(length(alltracks_mouse),2);
for i=1:length(trackingRes)
    if length(trackingFrameStamp{i})>=1000
        position(trackingFrameStamp{i},:)=trackingRes{i};
    end
end
position=smoothdata(position,'SmoothingFactor',0.01);   
position=position*behav.trackLength/behav.ROI3;
