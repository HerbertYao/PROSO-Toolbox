function [FBAsol,model_fit] = overlayTranscriptome(model,transcriptome,thres,waiver,method)

% A function that find a solution which fits a transcriptome count matrix
% as closely as possible. Reaction constraints have to be set up prior to
% using this function (such as minimal growth rate or production)
% 
% USAGE:
%   FBAsol = overlayTranscriptome(model_pc,count_matrix,5,[0.1,1000]);
% 
% INPUTS:
%   model:         A PC-model produced by function pcModel.m or preferably
%                  refined by adjustStoichAndKeff.m
%   transcriptome: Transcriptome data in count matrix form. A N*2 cell with
%                  geneId in column 1 and count in column 2
%   thres:         Lower threshold for a transcript count to be considered
%                  relevant
%   waiver:        A N*1 cell contains geneId that needs to be waivered
%                  from fitting
%   method:        'sqr' or 'abs'
%                       sqr: minimize the sum of square error
%                       abs: minimize the sum of absolute error
% 
% OUTPUTS:
%   FBAsol: FBA solution structure of best possible transcriptome fitting

if strcmpi(method,'abs')
%% Minimizing the Sum of Absolute Errors
    
%   Prepare model
    proteinExIdx = find(contains(model.rxns,'EX_protein_'));
    proteinList = model.mets(find(contains(model.mets,'protein_')));

%   Reset objective
    for i = 1:length(model.rxns)
        model.c(i) = 0;
    end

%   Add mets 'proteinFitHigh_' and 'proteinFitLow_' for each protein. 
%   It will be added to EX_protein reaction so generated together with 
%   protein conc.

    for i = 1:length(proteinExIdx)
        met = ['proteinFitHigh_',erase(proteinList{i},'protein_')];
        model = addMetabolite(model,met,'csense','L');
        model.S(length(model.mets),proteinExIdx(i)) = -1;

        met = ['proteinFitLow_',erase(proteinList{i},'protein_')];
        model = addMetabolite(model,met,'csense','G');
        model.S(length(model.mets),proteinExIdx(i)) = -1;

%       Changing proteinFit's b to count
%       proteinId string must either match or is contained in transcriptome
%       column 1. Otherwise b will be left 0
        idx = find(strcmp(transcriptome(:,1),erase(proteinList{i},'protein_')));
        if isempty(idx)
            idx = find(contains(transcriptome(:,1),erase(proteinList{i},'protein_')));
        end

        if length(idx) == 1
            if transcriptome{idx,2} >= thres
                model.b(length(model.mets)-1) = transcriptome{idx,2};
                model.b(length(model.mets)) = transcriptome{idx,2};
            end

        elseif isempty(idx) % if no count is found, proceed with warning
            warning('Count not found in transcriptome provided: %s\n',proteinList{i});

        else % if more than 1 count is found, proceed with warning
            if transcriptome{idx(1),2} >= thres
                model.b(length(model.mets)-1) = transcriptome{idx(1),2};
                model.b(length(model.mets)) = transcriptome{idx(1),2};
            end
            warning('Duplicate counts found in transcriptome provided: %s\n',proteinList{i});
        end

    end

%   Record their index
    proteinFitHighIdx = find(contains(model.mets,'proteinFitHigh_'));
    proteinFitLowIdx = find(contains(model.mets,'proteinFitLow_'));

%   Add proteinError_ reactions. They are the objective and need to be
%   minimized (min abs error)
    
    for i = 1:length(proteinFitHighIdx)
        model = addReaction(model,['proteinError_',erase(proteinList{i},'protein_')],...
            'metaboliteList',{model.mets{proteinFitHighIdx(i)},model.mets{proteinFitLowIdx(i)}},...
            'stoichCoeffList',[-1,1],...
            'reversible',false,...
            'objectiveCoef',1);
    end
    
%   Rescale proteinFitHigh_ and proteinFitLow_ .b so it sums to proteinWC.b
    sm = 0;
    for i = 1:length(proteinFitHighIdx)
        sm = sm + model.b(proteinFitHighIdx(i));
    end
    
    sc = model.b(find(strcmp(model.mets,'proteinWC')))/sm;
    
    for i = 1:length(proteinFitHighIdx)
        model.b(proteinFitHighIdx(i)) = model.b(proteinFitHighIdx(i)) * sc;
        model.b(proteinFitLowIdx(i)) = model.b(proteinFitLowIdx(i)) * sc;
    end
    
%   Giving waivers to waiverlist so their errors are not minimized
    if ~isempty(waiver)
        for i = 1:length(waiver)
            idx = find(strcmp(model.rxns,['proteinError_',waiver{i}]));
            model.c(idx) = 0;
        end
    end
    
%   Solve and return
    FBAsol = optimizeCbModel(model,'min');
    model_fit = model;

elseif strcmpi(method,'sqr')
%% Minimizing the Sum of Square Errors
    
%   initialize QP problem
    qp.A = model.S;
    qp.b = model.b;
    qp.c = zeros(length(model.c),1);
    qp.lb = model.lb;
    qp.ub = model.ub;
    qp.osense = 1;
    qp.csense = model.csense;
    qp.F = zeros(length(model.c),length(model.c));
    
%   Set qp.c and qp.F
    proteinExIdx = find(contains(model.rxns,'EX_protein_'));
    
    for i = 1:length(proteinExIdx)
        protName = erase(model.rxns{proteinExIdx(i)},'EX_');
        
%       Waiver
        if any(strcmp(waiver,protName))
            continue;
        end
        
%       Finding transcriptome entry and set them as qp.c
        idx = find(strcmp(transcriptome(:,1),erase(protName,'protein_')));
        if isempty(idx)
            idx = find(contains(transcriptome(:,1),erase(protName,'protein_')));
        end
        
        if length(idx) == 1 % one and only entry found in transcriptome
            if transcriptome{idx,2} >= thres
                qp.c(proteinExIdx(i)) = transcriptome{idx,2};
            end
            
        elseif isempty(idx) % no matched entry found
            warning('Count not found in transcriptome provided: %s\n',protName);
            
        else % multiple entries found
            if transcriptome{idx(1),2} >= thres
                qp.c(proteinExIdx(i)) = transcriptome{idx(1),2};
            end
            warning('Duplicate counts found in transcriptome provided: %s\n',protName);
        end
            
%       Also set qp.F. It is all 1 except proteins in the waiver list
        qp.F(proteinExIdx(i),proteinExIdx(i)) = 1;
    end
    
%   Rescale qp.c so it sums to protein budget
    sc = model.b(find(strcmp(model.mets,'proteinWC')))/sum(qp.c);
    for i = 1:length(qp.c)
        qp.c(i) = qp.c(i) * sc;
    end
    
%   Solve and prepare the return
    FBAsol = solveCobraQP(qp);
    model_fit = qp;
    
else
%% Error message
    error('Unrecognized method. Enter abs for absolute or sqr for square');
end

end