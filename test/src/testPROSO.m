% initCobraToolbox();
% initializePROSO();
% clc

inputs = struct();
outputs = struct();

% Test PC-FBA
% model: http://bigg.ucsd.edu/models/e_coli_core
fprintf('========== PC-FBA ==========\n');

model_ori = readCbModel('e_coli_core.xml');

inputs.model_ori = model_ori;

[model_pc_ori,fullProtein] = pcModel(model_ori, 200, 50);

FBAsol_pc = optimizeCbModel(model_pc_ori,'max');

outputs.FBAsol_pc = FBAsol_pc;

% Prepare RNA-seq data
% data: https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-GEOD-74882
fprintf('========== Data ==========\n');

expNames = ["E-GEOD-74882/GSM1937017_EV1_transcripts.txt",...
    "E-GEOD-74882/GSM1937018_EV2_transcripts.txt",...
    "E-GEOD-74882/GSM1937019_EV3_transcripts.txt",...
    "E-GEOD-74882/GSM1937020_EV4_transcripts.txt",...
    "E-GEOD-74882/GSM1937021_EV5_transcripts.txt",...
    "E-GEOD-74882/GSM1937022_EV6_transcripts.txt"];

for i = 1:length(expNames)
    d = readtable(expNames(i));
    [~, idx] = unique(d.Synonym);
    d = d(idx, :);
    if i == 1
        dt = d(:, [7,12]);
    else
        dt = join(dt, d(:, [7,12]));
    end
end
dt = dt(:, 2:end);

data = zeros(length(fullProtein),length(expNames));
waiver = {};

for i = 1:length(fullProtein)
    idx = find(contains(d.Synonym, fullProtein{i}));

    if ~isempty(idx)
        if length(idx) ~= 1
            warning('Protein %s: duplicating data entry\n',fullProtein{i});
        end
        data(i,:) = dt{idx,:};
    else
        waiver{end+1,1} = ['protein_',fullProtein{i}];
    end
end

% Test Non-convex QP
fprintf('========== Non-convex QP ==========\n');

inputs.model_pc_ori = model_pc_ori;
inputs.data = data;
inputs.waiver = waiver;

if ~exist('QPsol','var')
    [QPsol,model_qp] = overlayMultiomicsData(model_pc_ori,data,0,waiver,'keffEstimate',true);
    rIdx = find(startsWith(model_qp.varnames,'R_'));
    rValues = QPsol.x(rIdx);
end

outputs.rValues = rValues;

% Test Convex QP
fprintf('========== Convex QP ==========\n');

changeCobraSolverParams('QP','feasTol',1e-9);
changeCobraSolverParams('QP','printLevel',0);

[QPsol_cv,model_qp] = overlayMultiomicsData(model_pc_ori,data(:,1),0,waiver);
model_pc = implementProteinConstraints(model_pc_ori,QPsol_cv.full,waiver,0.02);
model_pc = model_pc{1};

outputs.QPsol_cv = QPsol_cv;

% Test Debottlenecker
fprintf('========== Debottleneck ==========\n');

inputs.model_pc = model_pc;

[FBAsol_db,model_db,relaxProt,relaxLevel] = proteinDebottleneck(model_pc,35);

outputs.FBAsol_db = FBAsol_db;

% Test Context-Specific FVA
fprintf('========== Context-Specific FVA ==========\n');

inputs.model_db = model_db;

FVAsol = contextSpecificPCFVA(model_db,inputs.model_ori.rxns,[0,0.5,0.9,0.99]);

outputs.FVAsol = FVAsol;

% Test PC-OptKnock
fprintf('========== PC-OptKnock ==========\n');

inputs.ok_target = 46;

[model_ok,param] = formulatePCOptKnock(model_pc_ori,1,0.05);
model_ok.obj(inputs.ok_target) = 1;
model_ok = changeNumKO(model_ok,3);
OKsol = gurobi(model_ok,param);

outputs.OKsol = OKsol;

% Test MOPA

inputs.ko_proteins = {'EX_protein_b1852', 'EX_protein_b2463', 'EX_protein_b1602'};

[MOPAsol,~,~] = MOPA(model_pc_ori,inputs.ko_proteins);

outputs.MOPAsol = MOPAsol;

% SAVE
save('../testData_testOfPROSO', 'inputs', 'outputs');
