function models = implementProteinConstraints(model,QPsols,waiverList,s)

if ~exist('s','var') % default of slack param
    s = 0.02;
end

proteinExIdx = find(startsWith(model.rxns,'EX_protein_'));
[~,w] = size(QPsols);

models = cell(w,1);

for i = 1:w

    model_alt = model; % create an unique instance for each experiment
    model_alt.lb(proteinExIdx) = QPsols(proteinExIdx,i);
    model_alt.ub(proteinExIdx) = QPsols(proteinExIdx,i) * (1-s);

    for j = 1:length(waiverList) % relief bounds for waiver list
        idx = find(strcmp(model.rxns,['EX_',waiverList{j}]));
        model_alt.lb(idx) = model.lb(idx);
        model_alt.ub(idx) = model.ub(idx);
    end

    models{i} = model_alt;
end

end