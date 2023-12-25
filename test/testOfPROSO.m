% The COBRAToolbox: testOfPROSO.m
%
% Purpose:
%     - Testing if all modules in the PROSO Toolbox is running properly
%
% Authors:
%     - Herbert Yao, Dec 16, 2023
%

global CBTDIR

% define the features required to run the test
requiredToolboxes = { ...
    'Statistics and Machine Learning Toolbox', ...
    'Bioinformatics Toolbox', ...
    'Deep Learning Toolbox' };

requiredSolvers = { 'gurobi' };

% require the specified toolboxes and solvers, along with a UNIX OS
% solversPkgs = prepareTest('requiredSolvers', requiredSolvers, 'requiredToolboxes', requiredToolboxes);
solversPkgs = prepareTest('requiredSolvers', requiredSolvers);

% save the current path and initialize the test
currentDir = cd(fileparts(which(mfilename)));

% determine the test path for references
testPath = pwd;

% set the tolerance
tol = 1e-4;

% load the model
model_ori = readCbModel('e_coli_core.xml');

%Load reference data
load([testPath, '/testData_testOfPROSO.mat']);

for k = 1:length(solversPkgs.LP)
    fprintf(' -- Running testOfPROSO using the solver interface: %s ... ', solversPkgs.LP{k});

    solverLPOK = changeCobraSolver(solversPkgs.LP{k}, 'all', 0);

    if solverLPOK

        % PC-FBA
        [model_pc_ori,fullProtein] = pcModel(model_ori, 200, 50);
        FBAsol_pc = optimizeCbModel(model_pc_ori,'max');
        if abs(outputs.FBAsol_pc.f - FBAsol_pc.f) > tol
            warning('PC-FBA solution exceeds tol');
        end

        % Non-convex QP
        if strcmp(solversPkgs.LP{k},'gurobi') % only proceed for gurobi
            [QPsol,model_qp] = overlayMultiomicsData(inputs.model_pc_ori,inputs.data,0,inputs.waiver,'keffEstimate',true);
            rValues = QPsol.x(find(startsWith(model_qp.varnames,'R_')));
            if abs(outputs.rValues(1) - rValues(1)) > tol
                warning('Non-convex QP solution exceeds tol');
            end
        else
            warning('Non-convex QP not tested for solver: %s', solversPkgs.LP{k});
        end

        % Convex QP
        [QPsol_cv,model_qp] = overlayMultiomicsData(inputs.model_pc_ori,inputs.data(:,1),0,inputs.waiver);
        if abs(outputs.QPsol_cv.full(52) - QPsol_cv.full(52)) > tol
            warning('Convex QP solution exceeds tol');
        end

        % Debottleneck
        [FBAsol_db,~,~,~] = proteinDebottleneck(inputs.model_pc,35);
        if abs(outputs.FBAsol_db.v(52) - FBAsol_db.v(52)) > tol
            warning('Debottleneck solution exceeds tol');
        end

        % Context-Specific FVA
        FVAsol = contextSpecificPCFVA(inputs.model_db,inputs.model_ori.rxns,[0,0.5,0.9,0.99]);
        if abs(outputs.FVAsol(52,1,8) - FVAsol(52,1,8)) > tol
            warning('Context-Specific FVA solution exceeds tol');
        end

        % Test PC-OptKnock
        [model_ok,param] = formulatePCOptKnock(inputs.model_pc_ori,1,0.05);
        model_ok = changeNumKO(model_ok,3);
        model_ok.c(inputs.ok_target) = 1;
        OKsol = solveCobraMILP(model_ok,param);
        if abs(outputs.OKsol.obj - OKsol.obj) > tol
            warning('PC-OptKnock solution exceeds tol');
        end

        % Test MOPA
        [MOPAsol,~,~] = MOPA(inputs.model_pc_ori,inputs.ko_proteins);
        if abs(outputs.MOPAsol.full(52) - MOPAsol.full(52)) > tol
            warning('MOPA solution exceeds tol');
        end

    end
%     wrongInputs = {'FirstArgument',modelForSecondArgument};
%     verifyCobraFunctionError('testFile', 'inputs', wrongInputs);
    % output a success message
    fprintf('Done.\n');
end

% change the directory
cd(currentDir)
