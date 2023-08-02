function [model_new,openRxnList] = tryCloseRxns(model,rxnList,closeDir,mu_min)
% Close as many designated reactions as possible while maintaining growth
% This is especially useful for setting up a medium with messy exchange rxn
% bounds.

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
