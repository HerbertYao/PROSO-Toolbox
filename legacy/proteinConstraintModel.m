function [constraintedModel,fullProtein,fullEnzyme,C_matrix,K_matrix,fullProteinMM] = proteinConstraintModel(model,fasta,maxWeightFrac)

% (MAIN) A function that produce a proteomically constrainted m-model from 
% a common m-model and protein sequence data. 
% 
% USAGE:
%   new_model = proteinConstraintModel(model,'Paeruginosa.fasta',550);
% 
% INPUTS:
%   model: A functional COBRA model struct with the field 'genes'
%   fasta: The file name of a proteome fasta file containing geneID 
%          (same set of ID as model.genes) as header and respective protein 
%          sequence in one-letter symbol.
%          Alternatively, you can put an estimated average protein length 
%          (avg # of amino acid per protein) if you don't have access to
%          fasta. The modelling accuracy will be compromised.
% 
% OPTIONAL INPUTS:
%   maxWeightFrac: A double denotes the maximum weight fraction of total
%                  protein components, in mg/gDW. Default = 550mg/gDW
%   keff_refTable: A n*2 cell array with enzyme ID in the first column and
%                  keff value in the second column, unit in 1/h. Enzymes 
%                  not found in this table will be assigned keff = 234000
%                  (NOT YET FINISHED)
% 
% OUTPUTS:
%   constraintedModel: A m-model with proteomic constraint
%   fullProtein:       A full list of proteins added to model.mets
%   fullEnzyme:        A full list of enzymes added to model.mets
%   C_matrix:          A matrix which related to enzyme formation
%   K_matrix:          A matrix of enzyme's catalyzing coefficient (k_eff)
%   fullProteinMM:     Molar mass of each protein in fullProtein
% 
% REQUIREMENTS:
%   MATLAB with Bioinformatics and Deep Learning Toolbox installed
%   Configured CobraToolbox with at least one capable solver
% 
% Important Note:
%   New_model.mets will contain all proteins and enzymes (protein
%   complexes). Each protein and enzyme will have an exchange rxn and a
%   enzyme formation rxn. Fluxes of these reactions are not real flux, but
%   representing their concentration in nmol/gDW. Therefore, the
%   concentration of protein_i is -v_{EX_protein_i} and the concentration 
%   of enzyme_j is v_{EX_enzyme_j} = v_{enzymeForm_enzyme_j}

% -------------------------------------------------------------------------
% Step 0: Parse Inputs

if isa(fasta,'char')
    fastaFile = fastaread(fasta);
    fprintf('Fasta file read successfully\n');
    
elseif isa(fasta,'double')
    avgProteinLen = fasta;
    fprintf('Input protein length: %d\n',avgProteinLen);
    
else
    error('Second input is neither an existing filename nor avg length\n');
end

if ~exist('maxWeightFrac','var')
    maxWeightFrac = 550;
    fprintf('Using default maximum protein weight fraction of 0.55\n');
    
elseif maxWeightFrac > 1000
    warning('Maximum proteome weight fraction is over 100%.\n');
    
elseif maxWeightFrac < 100
    warning('Maximum proteome weight fraction is under 10%.\n');
    
end

% -------------------------------------------------------------------------
% Step 1: List all unique enzymes and construct k_eff matrix
fprintf('Constructing k_eff matrix...');

% 1.1 List unique enzymes and put in the form of x(655)x(663)x(659)
fullEnzyme = {};

for i = 1:length(model.rules)
    if ~isempty(model.rules{i})
        try
            enzymeList = split(parseGeneRule(model.rules{i}),';');
        catch
            error('Error in geneRule parsing: rxn %d\n',i);
        end
        
%       Add new enzymes to fullEnzyme
        for j = 1:length(enzymeList)
            if ~any(strcmp(fullEnzyme,enzymeList{j}))
                fullEnzyme{length(fullEnzyme)+1,1} = enzymeList{j};
            end
        end
    end
end

% 1.2 Construct k_matrix such that vj <= sum_i (keff_ij*f)*ei
fullRxn = model.rxns;
K_matrix = zeros(length(fullRxn),length(fullEnzyme));

for i = 1:length(model.rules)
    enzymeList = split(parseGeneRule(model.rules{i}),';');
    
    for j = 1:length(enzymeList)
        idx = find(strcmp(fullEnzyme,enzymeList{j}));
        K_matrix(i,idx) = 234000;
    end
end
fprintf('done\n')

% -------------------------------------------------------------------------
% Step 2: Parse the fullEnzyme to construct C matrix
fprintf('Parsing geneRules and constructing C matrix...');

% Initialize C_matrix
fullProtein = model.genes;
C_matrix = zeros(length(fullProtein),length(fullEnzyme));

% Parse and assign values
for i = 1:length(fullEnzyme)
    
%   Remove brackets and the first 'x'
    enzyme = erase(fullEnzyme{i},{'(',')'});
    enzyme(1) = [];
    
%   Assign entries in C_matrix
    enzymeComp = split(enzyme,'x');
    for j = 1:length(enzymeComp)
        C_matrix(str2double(enzymeComp{j}),i) = 1;
    end
end

% May perform some checkings...
fprintf('done\n');

% -------------------------------------------------------------------------
% Step 3: Add proteins and enzymes to model.mets
%         Add EX_protein, EX_enzyme, and enzymeForm to model.rxns

% 3.1 Add protein and EX_protein first
fprintf('Adding protein metabolites and protein exchange reactions...');

for i = 1:length(fullProtein)
    name = ['protein_',fullProtein{i}];
    if ~any(strcmp(model.mets,name))
        model = addMetabolite(model,name);
        model = addExchangeRxn(model,name,-1000000,0);
    else
        warning('Replicate proteins in model: %s\n',name);
        model = addMetabolite(model,[name,'_1']);
        model = addExchangeRxn(model,[name,'_1'],-1000000,0);
    end
end

fprintf('done\n');

% 3.2 Add enzyme and EX_enzyme
fprintf('Adding enzymes metabolites and enzyme exchange reactions...');

for i = 1:length(fullEnzyme)
    name = ['enzyme_',fullEnzyme{i}];
    model = addMetabolite(model,name);
    model = addExchangeRxn(model,name,0,1000000);
end

fprintf('done\n');

% 3.3 Add enzymeForm rxns using C matrix
fprintf('Adding enzyme formation reactions...')

for i = 1:length(fullEnzyme)
    metList = {};
    coefList = [];
    
    for j = 1:length(fullProtein)
        if C_matrix(j,i) ~= 0
            metList{length(metList)+1} = ['protein_',fullProtein{j}];
            coefList(length(coefList)+1) = - C_matrix(j,i);
        end
    end
    
    metList{length(metList)+1} = ['enzyme_',fullEnzyme{i}];
    coefList(length(coefList)+1) = 1;
    name = ['enzymeForm_' fullEnzyme{i}];
    
    model = addReaction(model,name,'metaboliteList',metList,...
        'stoichCoeffList',coefList,'reversible',false);
end

fprintf('done\n');

% -------------------------------------------------------------------------
% Step 4: Add proteomic constraints to the model using K_matrix
fprintf('Adding proteomic constraints to the model...');

f = 1/1000000; % conversion factor from nmol to mmol

for i = 1:length(fullRxn)
    
%   Pass if the respective K_matrix row is all 0
    if ~any(K_matrix(i,:) ~= 0)
        if ~isempty(model.rules{i})
            warning('K matrix construction is wrong for rxn %s\n',model.rxns{i});
        end
        continue;
        
%   Else implement constraints
    else
        if isempty(model.rules{i})
            warning('K matrix construction is wrong for rxn %s\n',model.rxns{i});
        end
        constraintList = {};
        constraintCoef = [];
        constraintID = ['enzymeConstraint_',fullRxn{i}];
        
%       Look for non-zero keff and add respective enzyme into constraintList
        for j = 1:length(fullEnzyme)
            if K_matrix(i,j) ~= 0
                constraintList{length(constraintList)+1} = ['EX_enzyme_',fullEnzyme{j}];
                constraintCoef(length(constraintCoef)+1) = - K_matrix(i,j)*f;
            end
        end
        
%       Lastly add the constrainted rxn
        constraintList{length(constraintList)+1} = fullRxn{i};
        constraintCoef(length(constraintCoef)+1) = 1;
        
%       Implement constraints for the reaction
        model = addCOBRAConstraints(model,constraintList,0,'c',...
            constraintCoef,'dsense','L','ConstraintID',constraintID);
        
%       Also for the reverse reaction
        constraintID = [constraintID,'_REV'];
        constraintCoef(length(constraintCoef)) = -1;
        model = addCOBRAConstraints(model,constraintList,0,'c',...
            constraintCoef,'dsense','L','ConstraintID',constraintID);
    end
end

fprintf('done\n');

% -------------------------------------------------------------------------
% Step 5: Add total proteome weight constraint to the model
fprintf('Adding protein weight constraint to the model...');

% 5.1 Calculate each protein's molar mass either using avg length or fasta
fullProteinMM = zeros(length(fullProtein),1);

% If fasta is available, find the sequence first and calc molar mass
if exist('fastaFile','var')
    
%   Calculate the avg protein length
    avgLen = 0;
    for i = 1:length(fastaFile)
        avgLen = avgLen + length(fastaFile(i).Sequence);
    end
    avgLen = avgLen / length(fastaFile);
    
%   Find and calculate each protein's molar mass
    for i = 1:length(fullProteinMM)
        fullProteinMM(i) = calcProteinMM(findProteinSeq(fullProtein{i},fastaFile,avgLen));
    end
    
% If avgProteinLen is input by user, calculate the avgMM and assign to 
% every entry
elseif exist('avgProteinLen','var')
    avgSeq(1:avgProteinLen) = 'Z';
    avgMM = calcProteinMM(avgSeq);
    
    for i = 1:length(fullProteinMM)
        fullProteinMM(i) = avgMM;
    end
    
%   Else return error message
else
    error('Neither avg protein length or fasta file is found\n');
end

% 5.2 Implement the protein weight constraint

% Determining the constrainted 'rxn' list
fullProteinEx = reshape(fullProtein,[1,length(fullProtein)]);
for i = 1:length(fullProteinEx)
    fullProteinEx{i} = ['EX_protein_' fullProteinEx{i}];
end

% Implement
constraintedModel = addCOBRAConstraints(model,fullProteinEx,maxWeightFrac,...
    'c',transpose(-fullProteinMM*f),'dsense','L','ConstraintID','TotalProteinWeight');

fprintf('done\n');

end
