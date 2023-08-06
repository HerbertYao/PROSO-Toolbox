function [dFBAsol,substrateProfile,biomassProfile,model_new] = pcDynamicFBA(model,timeInt,bmIdx,initBM,substrateList,substrateConc,rbScale,model_ino)

% Dynamic FBA with ribosomal PC-model. 
% Principals are the same as regular dynamicFBA, but with extra 
% consideration over the total protein weight as well as the maximum 
% ribosome allocation.
% 
% USAGE:
% 
%   dFBAsol = pcDynamicFBA(pc_model,[0:1/60:1],1,{'EX_glc__D_e'},10,10000);
% 
% INPUTS:
% 
%   model:          A PC-model produced by function pcModel.m or preferably
%                   refined by adjustStoichAndKeff.m. All required exchange
%                   bounds need to be setup properly (including those in
%                   substrateList)
%   timeInt:        [N*1] Double denotes each time point when a sample is
%                   taken, unit in hour. Must start from zero
%   bmIdx:          The index of biomass reaction meant to be used
%   initBM:         Initial biomass in gDW/L
%   substrateList:  [M*1] Species of important substrates added to the 
%                   reactor. Exchange reactions for substrates are listed
%   substrateConc:  [M*1] Substrate concentrations in mol/L, matching 
%                   substrateList
% 
% OPTIONAL INPUPTS:
%   
%   model_ino: model used to simulate protein allocation at t = 0. Must be
%              the same model as the first input, but may with different 
%              lb and ub (to similate start-up lag phase)
% 
% NOT YET IMPLEMENTED:
%   noRibosome:     Number of ribosomes per gram dry weight
%   Vmax:           Maximum translation speed for each ribosome
% 
% OUTPUTS:
% 
%   dFBAsol:          Model solution as a function of time
%                     [length(model.rxns) x length(timeInt)]
%   substrateProfile: Substrates concentration as a function of time
%                     [length(substrateList) x length(timeInt)]
%   biomassProfile:   Biomass as a function of time
%                     [1 x length(timeInt)]

%% Step 0: Checking inputs

if length(timeInt) <= 1
    error('Time interval must be a double array with more than 1 element');
elseif timeInt(1) ~= 0
    error('Time interval must start from 0');
end

if initBM <= 0
    error('Initial biomass must be a positive value');
end

if any(substrateConc < 0)
    error('Substrate concentrations must be non-negative');
end

sol = optimizeCbModel(model,'max');
if sol.stat ~= 1
    error('The input model must be feasible');
end

for i = 1:length(substrateList)
    if ~any(strcmp(model.rxns,substrateList{i}))
        error('Model does not contain reaction: %s',substrateList{i});
    end
end

if length(substrateList) ~= length(substrateConc)
    error('Length of substrate list and substrate concentrations must be the same');
end

% Record length of model.rxns
rxnLen = length(model.rxns);
    
%% Step 1: Proceeding the first iteration for initial protein level

% Parse lists of protein from model
exProteinIdx = find(contains(model.rxns,'EX_protein_'));
proteinIdx = find(contains(model.mets,'protein_'));

if length(exProteinIdx) ~= length(proteinIdx)
    error('Number of proteins must equal to number of protein exchange reactions in PC model');
end

% Calculate and Record the initial protein allocation
if ~exist('model_ino','var')
    FBAsol = optimizeCbModel(model,'max');
    proteinConc0 = -FBAsol.v(exProteinIdx);
else
    FBAsol = optimizeCbModel(model_ino,'max');
    proteinConc0 = -FBAsol.v(exProteinIdx);
end

%% Step 2: Adding ribosomal constraint to model

% Need only a universal ribosome metabolite
model = addMetabolite(model,'ribosome','metName','ribosome resouces',...
    'b',1,'csense','L');

% Adding two new sets of constraints, proteinLB_ and proteinUB_
for i = 1:length(proteinIdx)
    model = addMetabolite(model,['proteinLB_',erase(model.mets{proteinIdx(i)},'protein_')],...
        'b',proteinConc0(i),'csense','G');
end

for i = 1:length(proteinIdx)
    model = addMetabolite(model,['proteinUB_',erase(model.mets{proteinIdx(i)},'protein_')],...
        'b',proteinConc0(i),'csense','L');
end

% Record their index
proteinLBIdx = find(contains(model.mets,'proteinLB_'));
proteinUBIdx = find(contains(model.mets,'proteinUB_'));

% Adding proteinLB and proteinUB to respective EX_protein rxns
for i = 1:length(proteinLBIdx)
    model.S(proteinLBIdx(i),exProteinIdx(i)) = -1;
    model.S(proteinUBIdx(i),exProteinIdx(i)) = -1;
end

% Implementing ribosome allocation to proteinLB and proteinUB
% In the form of ribAlloc
% THIS PART IS STILL TO BE FINISHED...CURRENTLY USABLE BUT VERY COARSE

proteinWCIdx = find(strcmp(model.mets,'proteinWC')); % for calculating dp

for i = 1:length(proteinIdx)
    
    dp = -model.S(proteinWCIdx,exProteinIdx(i))/rbScale; % FOR NOW USING SAME SCALE AS PROTEIN MW
    model = addReaction(model,['ribAlloc_',erase(model.mets{proteinIdx(i)},'protein_')],...
        'metaboliteList',{model.mets{proteinLBIdx(i)},model.mets{proteinUBIdx(i)},'ribosome'},...
        'stoichCoeffList',[1,-1,dp],...
        'reversible',false);
end

% Model setup completed

%% Step 3: Conducting dynamicFBA

% Configurate outputs
substrateProfile = zeros(length(substrateList),length(timeInt));
dFBAsol = zeros(rxnLen,length(timeInt));
biomassProfile = zeros(1,length(timeInt));

substrateProfile(:,1) = substrateConc;
biomassProfile(1) = initBM;

% Find idx for all substrate exchange rxn
substrateIdx = zeros(length(substrateList),1);
for i = 1:length(substrateIdx)
    substrateIdx(i) = find(strcmp(model.rxns,substrateList{i}));
end

% Start Iteration
for i = 1:length(timeInt)-1
    
%   delta_t for this iteration
    dt = timeInt(i+1)-timeInt(i);
    
%   Check if each substrate's concentration can sustain the next iteration.
%       C_subs[mol/L] - v_max[mol/h/gDW] * dt[h] * mu[gDW/L] >= 0
%   If not, change the bound.
    for j = 1:length(substrateList)
        if (substrateProfile(j,i) + model.lb(substrateIdx(j))*dt*biomassProfile(i)/1000 < 0)
            model.lb(substrateIdx(j)) = - substrateProfile(j,i)*1000/dt/biomassProfile(i);
        end
    end
    
%   Run FBA
    FBAsol = optimizeCbModel(model,'max');
    
%   If model infeasible at certain point, break and finish immediately
%   Don't return error because it's not necessarily a mistake
    if FBAsol.stat ~= 1
        warning('Model becomes infeasible at t = %.3f hr\n',timeInt(i));
        break;
        
%   If feasible, update changes
    else
%       Record all flux
        dFBAsol(:,i) = FBAsol.v(1:rxnLen);
        
%       Update biomass
        biomassProfile(i+1) = biomassProfile(i)*exp(FBAsol.v(bmIdx)*dt);
        
%       Update concentration of each substrate
        for j = 1:length(substrateList)
            substrateProfile(j,i+1) = substrateProfile(j,i) + FBAsol.v(substrateIdx(j))*dt*biomassProfile(i)/1000;
            if abs(substrateProfile(j,i+1)) < 1e-6 % NEED TO FIGURE OUT THIS NUMBER
                substrateProfile(j,i+1) = 0;
            end
        end
        
%       Update all proteinUB_ and proteinLB_
        for j = 1:length(proteinLBIdx)
            model.b(proteinLBIdx(j)) = -FBAsol.v(exProteinIdx(j));
            model.b(proteinUBIdx(j)) = -FBAsol.v(exProteinIdx(j));
        end
    end
    
end

model_new = model; % for debugging

end