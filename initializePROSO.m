function initializePROSO()

% Initialization and set up PROSO Toolbox

fprintf('============================\n');
fprintf('   PROSO Toolbox v1.0.0   \n');
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

% Add src to path and check the integrity of the file structure
lastwarn('','');
pth = erase(which('initializePROSO.m'),'initializePROSO.m');

addpath(pth);
addpath([pth,'src/']);
addpath([pth,'src/PC-FBA/']);
addpath([pth,'src/MOPA/']);
addpath([pth,'src/OVERLAY/']);
addpath([pth,'src/PC-OptKnock/']);
addpath([pth,'src/PC-DynamicFBA/']);
addpath([pth,'src/Utilities/']);

if ~isempty(lastwarn())
    error('PROSO Toolbox integrity problem');
end

end
