function [v1,v2] = drawSolutionEnvelope(model,rxnId1,rxnId2,NSteps,drawFlag,rxnWgt1,rxnWgt2)

if ~exist('rxnWgt1','var')
    rxnWgt1 = 1;
end
if ~exist('rxnWgt2','var')
    rxnWgt2 = 1;
end

% Determine rxn 1 flux range (slightly slacked)

model = changeObjective(model,rxnId1,rxnWgt1);
FBAsol = optimizeCbModel(model,'max');
v1_max = FBAsol.f - abs(FBAsol.f)*0.0001;
FBAsol = optimizeCbModel(model,'min');
v1_min = FBAsol.f + abs(FBAsol.f)*0.0001;

% Set v1 uniformly within [v1_min, v1_max]
% Also prealloc v2 vector

v1 = [v1_min:(v1_max-v1_min)/(NSteps-1):v1_max,...
    v1_max:(v1_min-v1_max)/(NSteps-1):v1_min]';

v2 = zeros(NSteps*2,1);

model = changeObjective(model,rxnId2,rxnWgt2);

% FVA

for i = 1:NSteps
    
    model = changeRxnBounds(model,rxnId1,v1(i),'b');
    FBAsol = optimizeCbModel(model,'max');
    v2(i) = FBAsol.f;
    FBAsol = optimizeCbModel(model,'min');
    v2(length(v2)-i+1) = FBAsol.f;

end

if drawFlag
    figure;
    plot(v1,v2,'LineWidth',3);
    xlabel(rxnId1);
    ylabel(rxnId2);
end

end