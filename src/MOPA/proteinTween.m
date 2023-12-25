function FBAsols = proteinTween(model,FBAsol_init,objRxn,nSteps)

% proteinTween finds the path from a certain solution to optima in nSteps.
% This is accomplished by tightly constraint the model's proteome to the
% previous step's solution and allows only a certain 'budget' of total 
% protein adjustment from the previous step (ribosomal PC-FBA). The output
% is the optimal trace of the cell moving from FBAsol_init to the optimum. 
% 
% USAGE:
% 
%   FBAsols = proteinTween(model_pc,FBAsol_old,'BIOMASS_SC5_notrace',20);
% 
% INPUTS:
%   model:       A PC-model formulated by function pcModel.m
%   FBAsol_init: The solution representing the current cell proteomic +
%                metabolic state. This can be a result from OVERLAY, MOPA,
%                or anything else. 
%   objRxn:      The objective reaction that is assumed to optimize over
%                a period of time. 
%   nSteps:      The number of steps taken between FBAsol_init and the
%                optimum. 
% 
% OUTPUTS:
%   FBAsols: A matrix with every ribosomal PC-FBA solutions. 
% 
% .. AUTHOR: Herbert Yao, Dec 2023
% 

model_ori = model;
proteinExIdx = find(startsWith(model.rxns,'EX_protein_'));

% Need to locate protein vectors for (1) the initial state and (2) final
% state. The final state is the closest protein vector to the initial while
% still allowing optimal growth 

if isfield(FBAsol_init,'v')
    protInit = FBAsol_init.v(proteinExIdx);
elseif isfield(FBAsol_init,'full')
    protInit = FBAsol_init.full(proteinExIdx);
elseif isfield(FBAsol_init,'x')
    protInit = FBAsol_init.x(proteinExIdx);
else
    error('Wrong FBAsol_init format');
end

% Change the objective

model = changeObjective(model,objRxn);

% Forcing the optimum and use convex QP to find the final protein state

FBAsol = optimizeCbModel(model,'max');
model = changeRxnBounds(model,objRxn,FBAsol.f,'b');

QPsol = overlayMultiomicsData(model,protInit); % note it's best to do unweighted here
protFin = QPsol.full(proteinExIdx);

% Record the grand difference between protFin and protInit

protDiff = sum(abs(protFin - protInit));

% Now setup a ribosomal PC-model

model = changeObjective(model_ori,objRxn);
model_ribo = formulateRibosomalPCModel(model,protDiff/nSteps);
model_ribo = updateRiboPCModel(model_ribo,FBAsol_init);

% Simulate each step

FBAsols = zeros(length(model.rxns),nSteps);

for i = 1:nSteps

    FBAsol = optimizeCbModel(model_ribo,'max');
    FBAsols(:,i) = FBAsol.v(1:length(model.rxns));
    model_ribo = updateRiboPCModel(model_ribo,FBAsol);
    
end

end