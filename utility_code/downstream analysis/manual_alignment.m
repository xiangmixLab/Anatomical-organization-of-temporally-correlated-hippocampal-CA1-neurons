function manual_alignment(fname)

%% step1: load neuron contours
all_Cn={};
for i=1:length(fname) % neurons, cells that contain neuron data
    load(fname{i});
    for j=1:length(neurons)
        all_Cn{i,j}=neurons{j}.Cn;
    end
end

%% step2: manual check aligned neurons and generate cell_register

all_cell_register={};
for i=1:length(fname) % neurons, cells that contain neuron data
    all_move_dis={};
    figure;
    for j=2:size(all_Cn,2)
        subplot(121)
        imagesc(all_Cn{i,1});
        subplot(122)
        imagesc(all_Cn{i,j});
        
        subplot(121)
        [x1,y1]=getpts; % 3 pts
        hold on;
        for k=1:length(x1)
            plot(x1,y1,'.','color','r','MarkerSize',6);
        end
        subplot(122)
        [x2,y2]=getpts; % 3 pts
        hold on;
        for k=1:length(x2)
            plot(x2,y2,'.','color','r','MarkerSize',6);
        end
        all_move_dis{j}=sum([x1-x2,y1-y2].^2,2).^0.5;
        close;
    end
    
end
neuron=Sources2D;
%Construct Sources2D object representing neurons over full recording
%This is really slow, optimize later

A_per_session=zeros(size(neurons{1}.A,1),size(cell_register,1),length(neurons));
for i=1:size(cell_register,1)
    
    A=zeros(size(neurons{1}.A(:,1)));
    
    C={};
    S={};
    C_raw={};
    for k=1:size(cell_register,2)
        if ~iszero(cell_register(i,k))
            C{k}=neurons{k}.C(cell_register(i,k),:);
            try
              C_raw{k}=neurons{k}.C_raw(cell_register(i,k),:);
            end
            try
                S{k}=neurons{k}.S(cell_register(i,k),:);
            end
            A=A+neurons{k}.A(:,cell_register(i,k));
            A_per_session(:,i,k)=neurons{k}.A(:,cell_register(i,k));
        else
            C{k}=zeros(1,size(neurons{k}.C,2));
            try
             C_raw{k}=zeros(1,size(neurons{k}.C,2));
            end
            try
                S{k}=zeros(1,size(neurons{k}.C,2));
            end
        end
    end
    
    
    A=A/size(cell_register,2);
    neuron.A=horzcat(neuron.A,A);
    neuron.C=vertcat(neuron.C,horzcat(C{:}));
    try
        neuron.C_raw=vertcat(neuron.C_raw,horzcat(C_raw{:}));
    end
    try
        neuron.S=vertcat(neuron.S,horzcat(S{:}));
    end
    centroid=calculateCentroid(A,data_shape(1),data_shape(2));
    neuron.centroid=vertcat(neuron.centroid,centroid);
    decays=[];
    try
    for q=1:size(cell_register,2)
        decays=[decays,neurons{q}.P.kernel_pars(cell_register(i,q))];
    end
 
    neuron.P.kernel_pars(i,1)=mean(decays);
    end
end
try
    neuron.Cn=neurons{base}.Cn;
end
neuron.imageSize=neurons{base}.imageSize;
neuron.updateCentroid();
try
    neuron=calc_snr(neuron);
end
% neuron.probabilities=aligned_probabilities;

neuron.cell_register=cell_register;
try
    neuron.options=neurons{1}.options;
end

%Eliminate remaining neurons falling below chain_prob
neuron.A_per_session=A_per_session;
neuron.delete(neuron.probabilities<chain_prob);
