function dpc_res=density_peak_clustering(img)

dpc_res=img*0;
radius=20;


for i=radius+1:size(img,1)-radius
    for j=radius+1:size(img,2)-radius
        dpc_res_t=0;
        for u=-radius:radius
            for v=-radius:radius
                dis2cen=((u-i)^2+(v-j)^2)^0.5;
                dpc_res_t=dpc_res_t+img(i,j)*(1/((2*pi)^0.5*radius))*exp(-dis2cen^2/(2*radius^2));
            end
        end
        dpc_res(i,j)=dpc_res_t;
    end
end