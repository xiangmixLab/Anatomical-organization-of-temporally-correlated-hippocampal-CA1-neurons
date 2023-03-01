% peak firing rate of field, cell series
% for olm, we suppose we already selected the obj proximal cells that the
% field with peak rate will be inside the field with obj
function peak_rates=peak_rate_field(fr)

for i=1:length(fr)
    