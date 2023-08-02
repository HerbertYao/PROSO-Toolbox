function model_ribo = formulateRibosomalPCModel(model,RibosomeBudget)
% This function formulates a ribosomal PC-model.
% Ribosomal PC-model aims to quantitatively constrain of amount of 
% proteome adjustments. In ribosomal PC-FBA, protein exchange reactions are
% effectively fixed at a certain value, and it must spend ribosome budget
% to deviate. 
% 
% USAGE:
% 
%   model_ribo = formulateRibosomalPCModel(model_pc,10);
% 
% INPUTS:
% 
%   model:          A PC-model formulated by function pcModel.m
%   RibosomeBudget: The total amount of protein adjustments allowed from
%                   the existing solution, unit in nmol/gDW. 
% 
% OUTPUTS:
% 
%   model_ribo: The ribosomal PC-model struct. Note that this model need to
%               be set by updateRiboPCModel.m before solving. The addtional
%               and modified mets and rxns of ribosomal PC-model are shown 
%               below. 
% 

% Format for ribosomal PC-model:
% 
%   New Mets LB_protein_b0001 and UB_protein_b0001:
%       LB_protein_b0001.b >= 0
%       UB_protein_b0001.b <= 0 (these will be set by updateRiboPCModel.m)
% 
%   New Met riboBudget:
%       riboBudget.b <= RibosomeBudget
% 
%   Modified Rxn EX_protein_b0001: 
%       protein_b0001 + 0.01 proteinWC + LB_protein_b0001 + UB_protein_b0001  <- 
% 
%   New Rxn ribo_protein_b0001:
%       UB_protein_b0001  ->  LB_protein_b0001 + riboBudget
% 

model_ribo = model;
proteinExIdx = find(startsWith(model.rxns,'EX_protein'));
proteinIdx = find(startsWith(model.mets,'protein_'));

% Add Met riboBudget

model_ribo = addMetabolite(model_ribo,'riboBudget','b',RibosomeBudget,'csense','L');

% Add Mets LB_protein and UB_protein

for i = 1:length(proteinExIdx)
    model_ribo = addMetabolite(model_ribo,['LB_',model.mets{proteinIdx(i)}],...
        'b',0,'csense','G');
    model_ribo = addMetabolite(model_ribo,['UB_',model.mets{proteinIdx(i)}],...
        'b',0,'csense','L');
end

lbProteinIdx = find(startsWith(model_ribo.mets,'LB_protein_'));
ubProteinIdx = find(startsWith(model_ribo.mets,'UB_protein_'));

% Modify Rxn EX_protein

for i = 1:length(proteinExIdx)
    model_ribo.S(lbProteinIdx(i),proteinExIdx(i)) = -1;
    model_ribo.S(ubProteinIdx(i),proteinExIdx(i)) = -1;
end

% Add Rxn ribo_protein

for i = 1:length(proteinExIdx)
    model_ribo = addReaction(model_ribo,['ribo_',model.mets{proteinIdx(i)}],...
        'metaboliteList',{model_ribo.mets{ubProteinIdx(i)},model_ribo.mets{lbProteinIdx(i)},'riboBudget'},...
        'stoichCoeffList',[-1 1 1],...
        'reversible',false);
end

end
