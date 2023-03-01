function img1=imadjustRGB2017(img)
img=uint8(img);
img(:,:,1)=imadjust(squeeze(img(:,:,1)));
img(:,:,2)=imadjust(squeeze(img(:,:,2)));
img(:,:,3)=imadjust(squeeze(img(:,:,3)));
img1=img;