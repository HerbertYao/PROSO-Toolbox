function initializeOVERLAY()

% Initialization and set up OVERLAY Toolbox

fprintf('============================\n');
fprintf('   OVERLAY Toolbox v1.0.1   \n');
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
pth = erase(which('initializeOVERLAY.m'),'initializeOVERLAY.m');

addpath(pth);
addpath([pth,'src/']);
addpath([pth,'src/PC-FBA/']);
addpath([pth,'src/MOPA/']);
addpath([pth,'src/OVERLAY/']);
addpath([pth,'src/PC-OptKnock/']);
addpath([pth,'src/Utilities/']);

if ~isempty(lastwarn())
    error('OVERLAY Toolbox integrity problem');
end

end
