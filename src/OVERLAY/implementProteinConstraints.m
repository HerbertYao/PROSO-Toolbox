function models = implementProteinConstraints(model,QPsols,waiverList,s)

% This function implements QP solutions to a PC-model.
% It produces a N*1 structure contains context-specific PC-models
% 
% USAGE:
% 
%   models = implementProteinConstraints(model_pc,QPsols,{'EX_protein_s0001'},0.03);
% 
% INPUTS:
% 
%   model:      A PC model built from the function 'pcModel.m'
%   QPsols:     The QP solution matrix constructed from concatenating
%               N number of convex QP solutions together. 
%   waiverList: This waiverList should be the same as used for the QP. 
%   s:          Universal slack term for the protein constraint. Please
%               check the tutorial or the publication for more details.
%               Default: 0.02
% 
% OUTPUTS:
% 
%   models: A N*1 struct containing context-specific PC-models 
% 

% default of slack param
if ~exist('s','var')
    s = 0.02;
end

% retain a reference model
if length(model) == 1
    model_ref = model;
else
    model_ref = model{1};
end

proteinExIdx = find(startsWith(model_ref.rxns,'EX_protein_'));
[~,w] = size(QPsols);

models = cell(w,1);

for i = 1:w

    if length(model) == 1
        model_alt = model; % create an unique instance for each experiment
    elseif length(model) == w
        model_alt = model{i};
    else
        error('Incorrect dimension for input: model');
    end

    model_alt.lb(proteinExIdx) = QPsols(proteinExIdx,i);
    model_alt.ub(proteinExIdx) = QPsols(proteinExIdx,i) * (1-s);

    for j = 1:length(waiverList) % relief bounds for waiver list
        idx = find(strcmp(model_ref.rxns,['EX_',waiverList{j}]));

        if length(model) == 1
            model_alt.lb(idx) = model.lb(idx);
            model_alt.ub(idx) = model.ub(idx);
        else
            model_alt.lb(idx) = model{i}.lb(idx);
            model_alt.ub(idx) = model{i}.ub(idx);
        end
    end

    models{i} = model_alt;
end

end