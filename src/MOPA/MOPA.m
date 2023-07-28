function [MOPAsol,FBAsol_wt,FBAsol_mut,matchTable] = MOPA(model,KOGene)

% Adapted from minimization of metabolic ajustment (MOMA), minimization of
% proteomic adjustment (MOPA) can find the sub-optimal solution when a cell
% just have a set of gene mutated based on its original mode of metabolism.
% A solver capable of solving QP is required
% 
% USAGE:
% 
%   FBAsol_mopa = MOPA(model_pc,{'EX_protein_s0001','EX_protein_s0002'});
% 
% INPUTS:
% 
%   model:    A PC-model produced by function pcModel.m or preferably
%             refined by adjustStoichAndKeff.m
%   delProts: Protein exchanges to be closed as a result of mutation
% 
% OUTPUTS:
% 
%   model_ok:
%   matchTable:
% 

% Obj: min (0.5*p'Ip - p_ori*p)
%   SSE = 2*Obj + p_ori'*I*p_ori
% 
% No need to exclude deleted genes from the objective as its error will be
% constant, although need to exclude them from the SSE calculation

% Define the QP problem
qp.A = model.S;
qp.b = model.b;
qp.c = zeros(length(model.c),1);
qp.lb = model.lb;
qp.ub = model.ub;
qp.osense = 1; % Unknown problem here: error if set to -1 (while all .c < 0). This should be identical
qp.csense = model.csense;
qp.F = zeros(length(model.c),length(model.c));

% find protein alloc from wt and assign them to qp.c
FBAsol = optimizeCbModel(model,'max');
proteinExIdx = find(contains(model.rxns,'EX_protein_'));

for i = 1:length(proteinExIdx)
    qp.c(proteinExIdx(i)) = -FBAsol.v(proteinExIdx(i));
end

% also assign qp.F = I only for p
for i = 1:length(proteinExIdx)
    qp.F(proteinExIdx(i),proteinExIdx(i)) = 1;
end

% Delete proteins in delGene
for i = 1:length(KOGene)
    idx = find(strcmp(model.rxns,KOGene{i}));
    qp.lb(idx) = 0;
    qp.ub(idx) = 0;
end

% Solve and return match table
MOPAsol = solveCobraQP(qp);
matchTable = zeros(length(proteinExIdx),2);

for i = 1:length(proteinExIdx)
    matchTable(i,1) = FBAsol.v(proteinExIdx(i));
    matchTable(i,2) = MOPAsol.full(proteinExIdx(i));
end

% Retrive wildtype optimum
FBAsol_wt = FBAsol;


% Tetrive KO strain optimum
model_mut = model;

for i = 1:length(KOGene)
    model_mut = changeRxnBounds(model_mut,KOGene{i},0);
end

FBAsol_mut = optimizeCbModel(model_mut,'max');

end
