function suppl17_h5_gen(dpath)

for i=1:length(dpath)
    ddata=[dpath{i},'\','final_concatenate_stack.mat'];
    load(ddata);
    Y(isnan(Y(:)))=0;
    Y(Y==inf)=0;
    Y(Y==-inf)=0;
    [ph,fn]=fileparts(ddata);
    h5create([ph,'\',fn,'.h5'],'/mov',size(Y))
    h5write([ph,'\',fn,'.h5'],'/mov',Y);
end