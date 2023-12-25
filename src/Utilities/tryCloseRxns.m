function [model_new,openRxnList] = tryCloseRxns(model,rxnList,closeDir,mu_min)

% Close as many designated reactions as possible while sustaining growth
% 
% This is especially useful for setting up a medium with messy exchange rxn
% bounds. For example, the usage example below shows how to find a feasible
% M-model with as little uptake reactions as possible. 
% 
% USAGE:
% 
%   [model_new,openRxnList] = tryCloseRxns(model,rxnList,closeDir,mu_min)
% 
% INPUTS: 
%   model:    An M-model or PC-model produced by function pcModel.m
%   rxnList:  The list of rxn to be closed
%   closeDir: Which way is each rxn in rxnList is shutting down. For
%             example, -1 for rxn "EX_o2_e" means we are attempting to make
%             an anaerobic environment (EX_o2_e.lb = 0), and 1 for the same
%             rxn means we want to stop the cell from producing oxygen
%             (EX_o2_e.ub = 0)
%   mu_min:   The minimal objective flux
% 
% OUTPUTS:
%   model_new:   The model with closed rxns lb or ub changed to 0
%   openRxnList: The list of rxns from rxnList that is left open
% 
% EXAMPLE:
% 
%   allExRxns = model_ori.rxns(find(startsWith(model_ori.rxns,'EX_')));
%   closeDir = -1 * ones(1,length(allExRxns));
%   model_new = tryCloseRxns(model_ori,allExRxns,closeDir,0.2);
% 
% .. AUTHOR: - Herbert Yao, Dec 2023
% 

model_milp = model;

model_milp.vartype = char('C' * ones(length(model_milp.rxns),1));

mu_idx = find(model_milp.c);
model_milp.lb(mu_idx) = mu_min;
model_milp.c(mu_idx) = 0;

for i = 1:length(rxnList)

    rxnIdx = find(strcmp(model_milp.rxns,rxnList{i}));

    model_milp = addMetabolite(model_milp,['BiCons_',rxnList{i}],...
        'b',0,'csense','G');
    model_milp.S(length(model_milp.mets),rxnIdx) = - closeDir(i);

    model_milp = addReaction(model_milp,['BiVar_',rxnList{i}],...
        'metaboliteList',{['BiCons_',rxnList{i}]},...
        'stoichCoeffList',1000,...
        'lowerBound',0,...
        'upperBound',1,...
        'objectiveCoef',1);

    model_milp.vartype(length(model.rxns)+i) = 'B';

end

model_milp.A = model_milp.S;
model_milp.osense = 1;

MILPsol = solveCobraMILP(model_milp);

openRxnList = rxnList(find(MILPsol.int));
model_new = model;

for i = 1:length(rxnList)
    if any(find(MILPsol.int) == i)
        continue;
    end

    if closeDir(i) == 1
        model_new.ub(find(strcmp(model_milp.rxns,rxnList{i}))) = 0;
    else
        model_new.lb(find(strcmp(model_milp.rxns,rxnList{i}))) = 0;
    end
end

end
