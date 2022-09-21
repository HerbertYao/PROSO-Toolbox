% ProteinConstraintModel testing

if ~exist('model_ori','var')
    model_ori = readCbModel('iGR774.xml');
end
    % Loading bounds data
if ~exist('boundData','var')
    addpath('pone.0229408.s002');
    [~, ~, boundData] = xlsread('pone.0229408.s005.xlsx');
end

model = model_ori;
% Setting up model
for i = 1:length(model_ori.rxns)
    
%   Looking up the reaction in the boundData
    idx = findBoundData(model.rxns{i},boundData);
    
%   First shutting down Nan and Pha specific rxns
    if (contains(model.rxns{i},'@Nan') && ~contains(model.rxns{i},'@Chl'))
        model = changeRxnBounds(model,model.rxns{i},0,'b');
        
    elseif (contains(model.rxns{i},'@Pha') && ~contains(model.rxns{i},'@Chl'))
        model = changeRxnBounds(model,model.rxns{i},0,'b');
        
%   Open Chl specific rxns and checking their reversibility
    elseif contains(model.rxns{i},'@Chl')
        model = changeRxnBounds(model,model.rxns{i},1000,'u');
        
%       Open lower bound if (1) rxn is found in table, and (2) is reversible
        if idx
            if str2double(boundData{idx,3})
                model = changeRxnBounds(model,model.rxns{i},-1000,'l');
            end
        end
        
%   Else check if the rxn appears in the boundData
    elseif idx
        model = changeRxnBounds(model,model.rxns{i},str2double(boundData{idx,4}),'l');
        model = changeRxnBounds(model,model.rxns{i},str2double(boundData{idx,5}),'u');
        
%   Else don't touch it for now
    else
%       model = changeRxnBounds(model,model.rxns{i},1000,'u');
%       model = changeRxnBounds(model,model.rxns{i},-1000,'l');
    end
end
clear idx
model = changeObjective(model,'@Chl_bio');

% Apply protein constraint
% model_pc = proteinConstraintModel(model,'all_fasta.fasta',600);

FBAsol = optimizeCbModel(model_pc,'max');

disp(FBAsol.f);


% =========================================================================
function idx = findBoundData(rxn, boundData)
    idx = 0;
    
    for i = 1:length(boundData)
        if strcmp(boundData{i,1},rxn)
            idx = i;
        end
    end
    
end
