%fix: Ntg or AD, or mice/cell idx
%randomm: randomm index, randomm, mice/cell index (can be multiple, in a
%cell)
%neuron: cell index

% for examples of formula, check matlab's fitlme document
function [lme,anova_lme]=calculate_lme_011221(var,fix,randomm,formulaa)

tbl_cell={};
vnames={'y'};

for i=1:length(fix)
    for j=1:length(var)
        tbl_cell{j,1}=var(j);
        tbl_cell{j,i+1}=fix{i}(j);
        vnames{i+1}=['fix',num2str(i)]; % fix1: first fix effect; fix2: second fix effect ...
    end
end

for i=1:length(randomm)
    for j=1:length(var)
        tbl_cell{j,i+length(fix)+1}=randomm{i}(j);
        vnames{i+length(fix)+1}=['random',num2str(i)]; % random1: first fix effect; random2: second fix effect ...
    end
end

tbl= cell2table(tbl_cell,'VariableNames',vnames);

lme = fitlme(tbl,formulaa); 

anova_lme=anova(lme);