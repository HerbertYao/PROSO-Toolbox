function initializePROSO(testPROSO)

% Initialization and set up PROSO Toolbox

fprintf('============================\n');
fprintf('    PROSO Toolbox v1.0.1    \n');
fprintf('============================\n');

% Check availability of COBRA Toolbox

if isempty(which('optimizeCbModel.m'))
    try
        warning('Boosting COBRA Toolbox\n');
        initCobraToolbox(false);
    catch
        error('COBRA Toolbox not installed or not in path');
    end
end

% Add src to path

PROSO_PATH = erase(which('initializePROSO.m'),'initializePROSO.m');

addpath(PROSO_PATH);
addpath(fullfile(PROSO_PATH,'src'));
addpath(fullfile(PROSO_PATH,'src','PC-FBA'));
addpath(fullfile(PROSO_PATH,'src','MOPA'));
addpath(fullfile(PROSO_PATH,'src','OVERLAY'));
addpath(fullfile(PROSO_PATH,'src','PC-OptKnock'));
addpath(fullfile(PROSO_PATH,'src','PC-DynamicFBA'));
addpath(fullfile(PROSO_PATH,'src','Utilities'));

% Check integrity of src

files = dir(fullfile(PROSO_PATH, '**', '*.m'));
fileNames = {files.name};

expectedFunctions = {
    'MOPA.m',... % MOPA
    'proteinTween.m',...
    'contextSpecificPCFVA.m',... % OVERLAY
    'implementProteinConstraints.m',...
    'overlayMultiomicsData.m',...
    'plotPCFVASolutions.m',...
    'proteinDebottleneck.m',...
    'updatePCModelKeffByR.m',...
    'formulateRibosomalPCModel.m',... % PC-DynamicFBA
    'pcDynamicFBA.m',...
    'updateRiboPCModel.m',...
    'adjustStoichAndKeff.m',... % PC-FBA
    'calcProteinMM.m',...
    'estimateKeffFromMW.m',...
    'findProteinSeq.m',...
    'parseGeneRule.m',...
    'pcModel.m',...
    'changeNumKO.m',... % PC-OptKnock
    'formulatePCOptKnock.m',...
};

missingFunctions = {};

for i = 1:numel(expectedFunctions)
    functionName = expectedFunctions{i};
    if ~any(strcmp(functionName, fileNames))
        missingFunctions = [missingFunctions, functionName];
    end
end

if ~isempty(missingFunctions)
    warning('Functions are missing from PROSO');
    fprintf('Missing functions: ');
    for i = 1:length(missingFunctions)
        fprintf('%s ',missingFunctions{i});
    end
    fprintf('\n');
end

% Run basic tests

if exist('testPROSO','var')
    if testPROSO == true
        currentDir = pwd;
        cd(fullfile(PROSO_PATH,'test'));
        testOfPROSO();
        cd(currentDir)
    end
end

end
