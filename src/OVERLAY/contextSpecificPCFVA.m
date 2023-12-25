function FVAsols = contextSpecificPCFVA(models_db,rxnList,optPerc,printBool)

% Conduct FVA for context-specific PC-models
% 
% USAGE:
%   
%   FVAsols = contextSpecificPCFVA(models_db,rxnList)
%   FVAsols = contextSpecificPCFVA(models_db,rxnList,optPerc)
%   FVAsols = contextSpecificPCFVA(models_db,rxnList,optPerc,printBool)
% 
% INPUTS:
%   models_db: M*1 cell contains context-specific PC-model
%   rxnList:   N*1 cell array with rxn IDs for FVA
% 
% OPTIONAL INPUTS:
%   optPerc: K*1 vector for objective flux optimal percentages. 
%            Default: [0,0.5,0.9]
%   printBool: if the coarse progress is printed in console.
% 
% OUTPUTS:
%   FVAsols: N * M * 2K matrix as FVA results. 
% 
% .. AUTHOR: - Herbert Yao, Dec 2023
% 

if ~exist('optPerc','var')
    optPerc = [0,0.5,0.9];
end
if ~exist('printBool','var')
    printBool = true;
end
if ~isa(models_db,'cell') && isa(models_db,'struct')
    models_db = {models_db};
end

FVAsols = zeros(length(rxnList),length(models_db),length(optPerc)*2);
if printBool
    fprintf('PC-FVA...');
end

for i = 1:length(models_db)
    if printBool
        fprintf('model %d/%d...',i,length(models_db));
    end
    model_alt = models_db{i};
    FBAsol = optimizeCbModel(model_alt,'max');
    objIdx = find(model_alt.c);

    for j = 1:length(optPerc)
        model_alt.lb(objIdx) = FBAsol.f * optPerc(j);

        for k = 1:length(rxnList)
            model_alt = changeObjective(model_alt,rxnList{k});
            FBAsol_alt = optimizeCbModel(model_alt,'min');
            FVAsols(k,i,j) = FBAsol_alt.f;
        end
    end

    for j = 1:length(optPerc)
        model_alt.lb(objIdx) = FBAsol.f * optPerc(length(optPerc)+1-j);

        for k = 1:length(rxnList)
            model_alt = changeObjective(model_alt,rxnList{k});
            FBAsol_alt = optimizeCbModel(model_alt,'max');
            FVAsols(k,i,j+4) = FBAsol_alt.f;
        end
    end
end

if printBool
    fprintf('done\n');
end

end