function [pcOptKnockMILP,param] = formulatePCOptKnock(model,numKO,objMin,FVABool,PCBool)
%%
% This function prepares for pcOptKnock by formulating a MILP to be fed in
% the optKnock
% 
% USAGE:
% 
%   milp_ok = formulatePCOptKnock(model_pc);
% 
% INPUTS:
% 
%   model:       A PC-model produced by function pcModel.m or preferably
%                refined by adjustStoichAndKeff.m
%   numKO: The maximum number of knocked out allowed. This can be later
%   changed by changeNumKO
%   
% 
% OUTPUTS:
% 
%   pcOptKnockMILP: MILP struct formulated for Gurobi optimizer
%   param:          Solver parameters that can be passed in to the solver
% 

%% Prepare original MILP from model
if ~exist('numKO','var')
    numKO = 1;
end
if ~exist('objMin','var')
    objMin = 0.05;
end
if ~exist('FVABool','var')
    FVABool = true;
end
if ~exist('PCBool','var')
    PCBool = true;
end

fprintf('Configuring MILP from model...');

% Further constrain PC-model
if FVABool
    fprintf('FVA...');
    for i = 1:length(model.rxns)
        model_alt = changeObjective(model,model.rxns{i});
        FBAsol = optimizeCbModel(model_alt,'max');
        model.ub(i) = FBAsol.f;
        FBAsol = optimizeCbModel(model_alt,'min');
        model.lb(i) = FBAsol.f;
    end
else
    for i = 1:length(model.rxns)
        if model.lb(i) < -1000
            model.lb(i) = -1000;
        end
        if model.ub(i) > 1000
            model.ub(i) = 1000;
        end
    end
end

% Modify proteinWC so its b is also zero. Also all csense are now E
if PCBool
    model = addExchangeRxn(model,'proteinWC',0,model.b(find(strcmp(model.mets,'proteinWC'))));
    model.b(find(strcmp(model.mets,'proteinWC'))) = 0;
    model.csense(find(strcmp(model.mets,'proteinWC'))) = 'E';
    model.csense(find(startsWith(model.mets,'protein_'))) = 'E';
end

% Define MILP
milp.A = model.S;
milp.rhs = model.b;
milp.obj = zeros(length(model.c),1);
milp.lb = model.lb;
milp.ub = model.ub;
milp.modelsense = 'max';
milp.vtype = char(ones(length(model.rxns),1)*'C');

milp.sense = model.csense;
milp.sense(milp.sense == 'E') = '=';
milp.sense(milp.sense == 'G') = '>';
milp.sense(milp.sense == 'L') = '<';

% Keep rxn names and met names for the sake of finding entries
milp.varnames = model.rxns;
milp.constrnames = model.mets;

fprintf('done\n');

%% Bilevel Variables
fprintf('Adding bilevel variables...');

% Add variable w for each metabolite
for i = 1:length(model.mets)
    [hei,wid] = size(milp.A);
    milp.A(:,wid+1) = zeros(hei,1);
    milp.obj(wid+1) = 0;
    milp.lb(wid+1) = -1000;
    milp.ub(wid+1) = 1000;
    milp.vtype(wid+1) = 'C';
    milp.varnames{wid+1} = ['w_',model.mets{i}];
end

% Add variable a, b, sigma, and y for each reaction
for i = 1:length(model.rxns)
%   a
    [hei,wid] = size(milp.A);
    milp.A(:,wid+1) = zeros(hei,1);
    milp.obj(wid+1) = 0;
    milp.lb(wid+1) = 0;
    milp.ub(wid+1) = 1000;
    milp.vtype(wid+1) = 'C';
    milp.varnames{wid+1} = ['a_',model.rxns{i}];
    
%   b
    [hei,wid] = size(milp.A);
    milp.A(:,wid+1) = zeros(hei,1);
    milp.obj(wid+1) = 0;
    milp.lb(wid+1) = 0;
    milp.ub(wid+1) = 1000;
    milp.vtype(wid+1) = 'C';
    milp.varnames{wid+1} = ['b_',model.rxns{i}];
    
%   y
    [hei,wid] = size(milp.A);
    milp.A(:,wid+1) = zeros(hei,1);
    milp.obj(wid+1) = 0;
    milp.lb(wid+1) = 0;
    milp.ub(wid+1) = 1;
    milp.vtype(wid+1) = 'I';
    milp.varnames{wid+1} = ['y_',model.rxns{i}];
    
%   sigma
    [hei,wid] = size(milp.A);
    milp.A(:,wid+1) = zeros(hei,1);
    milp.obj(wid+1) = 0;
    milp.lb(wid+1) = -1000;
    milp.ub(wid+1) = 1000;
    milp.vtype(wid+1) = 'C';
    milp.varnames{wid+1} = ['sigma_',model.rxns{i}];
end

% Record a list of all variables added
wIdx = find(startsWith(milp.varnames,'w_'));
aIdx = find(startsWith(milp.varnames,'a_'));
bIdx = find(startsWith(milp.varnames,'b_'));
yIdx = find(startsWith(milp.varnames,'y_'));
sigmaIdx = find(startsWith(milp.varnames,'sigma_'));

fprintf('done\n');

%% Bilevel Constraints
fprintf('Adding bilevel constraints...');

% Dual constraints for each reaction in the model
fprintf('dual cons...');
% =========================
%   S'w + b - a + sigma = c
% =========================

for i = 1:length(model.rxns)
    [hei,wid] = size(milp.A);
    
    milp.A(hei+1,:) = zeros(1,wid);
%   S'w
    idx = find(model.S(:,i));
    for j = 1:length(idx)
        milp.A(hei+1,wIdx(idx(j))) = model.S(idx(j),i);
    end
%   1*b
    milp.A(hei+1,bIdx(i)) = 1;
%   -1*a
    milp.A(hei+1,aIdx(i)) = -1;
%   1*sigma
    milp.A(hei+1,sigmaIdx(i)) = 1;
    
    milp.rhs(hei+1) = model.c(i); % = c
    milp.sense(hei+1) = '=';
    milp.constrnames{hei+1} = ['dual_',model.rxns{i}];
end

fprintf('mu/eta...');
% ================
% a - mu_U*y <= 0
% b - eta_U*y <= 0
% ================

for i = 1:length(model.rxns)
    [hei,wid] = size(milp.A);
    milp.A(hei+1,:) = zeros(1,wid);
    milp.A(hei+1,aIdx(i)) = 1;
    milp.A(hei+1,yIdx(i)) = -1000;
    milp.rhs(hei+1) = 0;
    milp.sense(hei+1) = '<';
    milp.constrnames{hei+1} = ['mu_U_',model.rxns{i}];
    
    [hei,wid] = size(milp.A);
    milp.A(hei+1,:) = zeros(1,wid);
    milp.A(hei+1,bIdx(i)) = 1;
    milp.A(hei+1,yIdx(i)) = -1000;
    milp.rhs(hei+1) = 0;
    milp.sense(hei+1) = '<';
    milp.constrnames{hei+1} = ['eta_U_',model.rxns{i}];
end

fprintf('sigma bounds...');
% ub and lb on sigma
% ================
% sigma - M*y >= -M
% sigma + M*y <= M
% ================

for i = 1:length(model.rxns)
    [hei,wid] = size(milp.A);
    milp.A(hei+1,:) = zeros(1,wid);
    milp.A(hei+1,sigmaIdx(i)) = 1;
    milp.A(hei+1,yIdx(i)) = -1000;
    milp.rhs(hei+1) = -1000;
    milp.sense(hei+1) = '>';
    milp.constrnames{hei+1} = ['sigma_L_',model.rxns{i}];
    
    [hei,wid] = size(milp.A);
    milp.A(hei+1,:) = zeros(1,wid);
    milp.A(hei+1,sigmaIdx(i)) = 1;
    milp.A(hei+1,yIdx(i)) = 1000;
    milp.rhs(hei+1) = 1000;
    milp.sense(hei+1) = '<';
    milp.constrnames{hei+1} = ['sigma_U_',model.rxns{i}];
end

fprintf('strongDuality...');
% set primal obj = dual obj for optimizing objective (strong duality)
% ===================================
% 0*w + u'b - l'a + 0*sigma - c'v = 0
% ===================================

[hei,wid] = size(milp.A);
milp.A(hei+1,:) = zeros(1,wid);

for i = 1:length(model.rxns)
    milp.A(hei+1,bIdx(i)) = model.ub(i);
    milp.A(hei+1,aIdx(i)) = -model.lb(i);
    milp.A(hei+1,i) = -model.c(i);
end

milp.rhs(hei+1) = 0;
milp.sense(hei+1) = '=';
milp.constrnames{hei+1} = 'strongDuality';

fprintf('binary cons...');
% ===========
% v - ly >= 0
% v - uy <= 0
% ===========

for i = 1:length(model.rxns)
    [hei,wid] = size(milp.A);
    milp.A(hei+1,:) = zeros(1,wid);
    milp.A(hei+1,i) = 1;
    milp.A(hei+1,yIdx(i)) = -model.lb(i);
    milp.rhs(hei+1) = 0;
    milp.sense(hei+1) = '>';
    milp.constrnames{hei+1} = ['lb_y_',model.rxns{i}];
    
    [hei,wid] = size(milp.A);
    milp.A(hei+1,:) = zeros(1,wid);
    milp.A(hei+1,i) = 1;
    milp.A(hei+1,yIdx(i)) = -model.ub(i);
    milp.rhs(hei+1) = 0;
    milp.sense(hei+1) = '<';
    milp.constrnames{hei+1} = ['ub_y_',model.rxns{i}];
end

fprintf('Max KO (K = %d)...',numKO);
% ===============
% sum(y) >= N - K
% ===============

[hei,wid] = size(milp.A);

milp.A(hei+1,:) = zeros(1,wid);
for i = 1:length(yIdx)
    milp.A(hei+1,yIdx(i)) = 1;
end

milp.rhs(hei+1) = length(yIdx) - numKO; % Allowing 1 KO by default
milp.sense(hei+1) = '>';
milp.constrnames{hei+1} = 'max_KO';

% ==========================
% v_biomass >= v_biomass^min
% ==========================

milp.lb(find(model.c)) = objMin;

fprintf('done\n');

%% Adjust binary variables

% Only allow y_protein to be knockout
if PCBool
    fprintf('Fixing non-protein binary variables to 1...');
    
    for i = 1:length(yIdx)
        if ~contains(milp.varnames{yIdx(i)},'y_EX_protein_')
            milp.lb(yIdx(i)) = 1;
        end
    end
end

% Fixing binary vars of growth-essential proteins to 1
if PCBool && FVABool
    fprintf('Keeping growth-essential proteins open...');
    
    for i = 1:length(model.rxns)
        if contains(model.rxns{i},'EX_protein_')
            
            model_alt = model;
            model_alt.lb(i) = 0;
            FBAsol_alt = optimizeCbModel(model_alt,'max');
            
            if (FBAsol_alt.stat ~= 1) || (FBAsol_alt.f < 1e-6)
                idx = find(strcmp(milp.varnames,['y_',model.rxns{i}]));
                milp.lb(idx) = 1;
            end
        end
    end
    
    fprintf('done\n');
end

if ~issparse(milp.A)
    milp.A = sparse(milp.A);
end

pcOptKnockMILP = milp;
param.FeasibilityTol = 1e-9;
param.IntFeasTol = 1e-9;

end