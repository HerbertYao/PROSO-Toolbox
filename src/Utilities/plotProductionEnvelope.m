function [v1,v2] = plotProductionEnvelope(models,rxnId1,rxnId2,NSteps,plotFlag,lgd,rxnWgt1,rxnWgt2)

% This function draws a 2d production envelope for each PC-model
% 
% USAGE:
% 
%   plotProductionEnvelope(model_pc,'BIOMASS_SC5_notrace','EX_succ_e',20,true);
%   plotProductionEnvelope(models,'BIOMASS_SC5_notrace','EX_succ_e',20,true);
% 
% INPUTS:
% 
%   models:   A PC-model formulated by function pcModel.m
%   rxnId1:   The first rxn of interests
%   rxnId2:   The second rxn of interests
%   NSteps:   Number of points taken between v1_max and v1_min
%   drawFlag: If the function plot a figure
%   Lgd:      Legend for the plot. Must has the same length as models
%   rxnWgt1:  % TODO %
%   rxnWgt2:  % TODO %
% 
% OUTPUTS:
% 
%   v1: Flux vector for rxnId1
%   v2: Flux vector for rxnId2
% 

if ~exist('rxnWgt1','var')
    rxnWgt1 = 1;
end
if ~exist('rxnWgt2','var')
    rxnWgt2 = 1;
end

% For each model in models

for j = 1:length(models)

    if length(models) > 1
        model = models{j};
    else
        model = models;
    end

    % Determine rxn 1 flux range (slightly slacked)
    
    model = changeObjective(model,rxnId1);
    FBAsol = optimizeCbModel(model,'max');
    v1_max = FBAsol.f - abs(FBAsol.f)*0.0001;
    FBAsol = optimizeCbModel(model,'min');
    v1_min = FBAsol.f + abs(FBAsol.f)*0.0001;
    
    % Set v1 uniformly within [v1_min, v1_max]
    % Also prealloc v2 vector
    
    if j == 1
        v1 = [v1_min:(v1_max-v1_min)/(NSteps-1):v1_max,...
            v1_max:(v1_min-v1_max)/(NSteps-1):v1_min]';
        v2 = zeros(NSteps*2,1);
    else
        v1 = [v1, [v1_min:(v1_max-v1_min)/(NSteps-1):v1_max,...
            v1_max:(v1_min-v1_max)/(NSteps-1):v1_min]'];
        v2 = [v2, zeros(NSteps*2,1)];
    end
    
    model = changeObjective(model,rxnId2);
    
    % FVA
    
    for i = 1:NSteps
        
        model = changeRxnBounds(model,rxnId1,v1(i,j),'b');
        FBAsol = optimizeCbModel(model,'max');
        v2(i,j) = FBAsol.f;
        FBAsol = optimizeCbModel(model,'min');
        v2(length(v2)-i+1,j) = FBAsol.f;
    
    end
end

% Plot 

if plotFlag
    figure;
    plot(v1,v2,'LineWidth',3);
    xlabel(rxnId1,'Interpreter','none');
    ylabel(rxnId2,'Interpreter','none');
    if exist('lgd','var')
        legend(lgd);
    end
end

end
