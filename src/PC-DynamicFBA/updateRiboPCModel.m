function model_ribo_new = updateRiboPCModel(model_ribo,FBAsol)

% This function update ribosomal PC-model's existing proteome vector. 
% 
% USAGE:
% 
%   model_ribo = updateRiboPCModel(model_ribo,FBAsol);
% 
% INPUTS:
%   model_ribo: A ribosomal PC-model from formulateRibosomalPCModel.m
%   FBAsol:     The LP, MILP, or QP solution representing the cell's 
%               previous operating point
% 
% OUTPUTS:
%   model_ribo_new: The updated ribosomal PC-model
%  
% .. AUTHOR: Herbert Yao, Dec 2023
% 

proteinExIdx = find(startsWith(model_ribo.rxns,'EX_protein_'));
lbProteinIdx = find(startsWith(model_ribo.mets,'LB_protein_'));
ubProteinIdx = find(startsWith(model_ribo.mets,'UB_protein_'));

if isfield(FBAsol,'v')
    model_ribo.b(lbProteinIdx) = -FBAsol.v(proteinExIdx)*0.999;
    model_ribo.b(ubProteinIdx) = -FBAsol.v(proteinExIdx);
elseif isfield(FBAsol,'x')
    model_ribo.b(lbProteinIdx) = -FBAsol.x(proteinExIdx)*0.999;
    model_ribo.b(ubProteinIdx) = -FBAsol.x(proteinExIdx);
elseif isfield(FBAsol,'full')
    model_ribo.b(lbProteinIdx) = -FBAsol.full(proteinExIdx)*0.999;
    model_ribo.b(ubProteinIdx) = -FBAsol.full(proteinExIdx);
else
    error('Wrong FBAsol format');
end

model_ribo_new = model_ribo;

end