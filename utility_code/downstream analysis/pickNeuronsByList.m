function neuron=pickNeuronsByList(neuron0,cellList)

try 
    neuron=neuron0.copy;
catch
    neuron=neuron0;
end

if isfield(neuron,'A')
    neuron.A=neuron.A(:,cellList);
end
if isfield(neuron,'C')
    neuron.C=neuron.C(cellList,:);
end
if isfield(neuron,'S')
    neuron.S=neuron.S(cellList,:);
end
if isfield(neuron,'trace')
    neuron.trace=neuron.trace(cellList,:);
end
if isfield(neuron,'C_raw')
    neuron.C_raw=neuron.C_raw(cellList,:);
end