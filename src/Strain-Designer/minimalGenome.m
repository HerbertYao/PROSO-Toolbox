function [model_min,KOList] = minimalGenome(model,percentOri,protectList)

% Finding a minimal cell by shutting down as many protein exchanges as
% possible while keep the objective flux at a certain percentage of the
% original
% 
% USAGE:
% 
%   model_min = findMinimalCell(model_pc,0.5);
% 
% INPUTS:
% 
%   model:       A PC-model produced by function pcModel.m or preferably
%                refined by adjustStoichAndKeff.m
%   percentOri:  
%   protectList: List of proteins mets that are protected from deletion
% 
% OUTPUTS:
% 
%   model_min:
%   KOList:
% 

model_min = model; % Keep the original model

% Put a constraint on the original objective, then remove its objective
% coefficient
FBAsol = optimizeCbModel(model,'max');
model.lb(find(model.c)) = FBAsol.f * percentOri;
model.c(find(model.c)) = 0;

% Add yu_protein metabolites to EX_protein reactions
proteinExIdx = find(contains(model.rxns,'EX_protein_'));

for i = 1:length(proteinExIdx)
    model = addMetabolite(model,['yu_',erase(model.rxns{proteinExIdx(i)},'EX_')],...
        'csense','G');
    model.S(length(model.mets),proteinExIdx(i)) = 1;
end

% Add binary variables which supplying yu_protein
% They are the new objective
for i = 1:length(proteinExIdx)
    model = addReaction(model,['bvar_',erase(model.rxns{proteinExIdx(i)},'EX_')],...
        'metaboliteList',{['yu_',erase(model.rxns{proteinExIdx(i)},'EX_')]},...
        'stoichCoeffList',[1000000],...
        'lowerBound',0,...
        'upperBound',1,...
        'objectiveCoef',1);
end

bvarIdx = find(contains(model.rxns,'bvar_'));

% Protect proteins in the protectList
if exist('protectList','var')
    for i = 1:length(protectList)
        idx = find(strcmp(model.rxns,['bvar_',protectList{i}]));
        model.lb(idx) = 1;
    end
end

% Config the MILP
milp.A = model.S;
milp.b = model.b;
milp.c = model.c;
milp.lb = model.lb;
milp.ub = model.ub;
milp.osense = 1; % Minimize the number of proteins openned
milp.csense = model.csense;
milp.x0 = zeros(length(model.c),1);
milp.vartype = model.csense;

% Specify binary variables
for i = 1:length(model.rxns)
    if any(bvarIdx == i)
        milp.vartype(i) = 'B';
    else
        milp.vartype(i) = 'C';
    end
end

% Run MILP solver and collect result
sol = solveCobraMILP(milp,'intTol',1e-8);

KOList = {};

for i = 1:length(bvarIdx)
    if sol.full(bvarIdx(i)) == 0
        model_min.lb(proteinExIdx(i)) = 0;
        KOList{length(KOList)+1,1} = erase(model_min.rxns{proteinExIdx(i)},'EX_');
    end
end

end