function dis_ang_plot=obj_dis_ang_plot_convert(fr,obj)

origin_p=[size(fr,2)/2,size(fr,1)/2]; % x,y

dis_list=[];
ang_list=[];
corresponding_fr_idx=[];
ctt=1;
for i=1:size(fr,1)
    for j=1:size(fr,2)
       dis_list(ctt)=round(norm([[j,i]-origin_p;obj-origin_p]'));
       ang_list(ctt)= round(vectorAngles(obj-origin_p,[j,i]-origin_p));
       corresponding_fr_idx(ctt)=sub2ind(size(fr),i,j);
       ctt=ctt+1;
    end
end

dis_ang_plot=[];
for i=1:length(dis_list)    
   dis_ang_plot(dis_list(i)+1,ang_list(i)+181)=fr(corresponding_fr_idx(i));
end

dis_ang_plot=filter2DMatrices(dis_ang_plot,1);