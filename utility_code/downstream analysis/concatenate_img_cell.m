function vid=concatenate_img_cell(V,toGraySign)

vid=[];
if length(size(V{1}))==2
    for i=1:length(V)
        vid(:,:,i)=V{i};
    end    
end

if length(size(V{1}))==3
    if toGraySign==0
        for i=1:length(V)
            vid(:,:,:,i)=V{i};
        end   
    else
        for i=1:length(V)
            vid(:,:,i)=rgb2gray(V{i});
        end  
    end
end