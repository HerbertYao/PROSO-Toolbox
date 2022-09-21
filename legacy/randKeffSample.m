function [result_keff,result_rxn,result_enzyme,result_protein] = randKeffSample(model,noSample,lb,ub)

% Function that samples all k_eff randomly for a protein-constrained
% m-model
% 
% USAGE:
%   [result_keff,result_v,result_enzyme,result_protein] = randKeffSample(model,500,108000,360000);
% 
% INPUTS:
%   model:    A protein-constrained COBRA model
%   noSample: No. of samples to be taken
%   lb:       Lower bound of randomized K_eff
%   ub:       Upper bound of randomized K_eff
% 
% OUTPUTS:
%   result_eff:     Matrix of 
%   result_v:       A full list of proteins added to model.mets
%   result_enzyme:  A full list of enzymes added to model.mets
%   result_protein: A matrix which related to enzyme formation
% 

% List real reactions, protein exchanges, and enzyme exchanges
f = 1e-6;
fullRxn = {};
fullProtein = {};
fullEnzyme = {};

for i = 1:length(model.rxns)
    if contains(model.rxns{i},'EX_protein_')
        fullProtein{length(fullProtein)+1,1} = model.rxns{i};
    elseif contains(model.rxns{i},'EX_enzyme_')
        fullEnzyme{length(fullEnzyme)+1,1} = model.rxns{i};
    elseif contains(model.rxns{i},'enzymeForm_')
        continue;
    else
        fullRxn{length(fullRxn)+1,1} = model.rxns{i};
    end
end

% Allocate outputs
result_keff = zeros(length(fullEnzyme),noSample);
result_rxn = zeros(length(fullRxn),noSample);
result_protein = zeros(length(fullProtein),noSample);
result_enzyme = zeros(length(fullEnzyme),noSample);

% Sampling
fprintf('\nSampling...');
for i = 1:noSample
%   Generate random k_eff values
    k_eff = random('Uniform', lb, ub, [length(fullEnzyme),1]);
    
%   Assign k_eff to PC model
%       model.C contains constraint coefficients
    for j = 1:(length(model.ctrs)-1) % Exclude the last constraint (which is the total enzyme mass)
        
%       Protein constraint was applied to EX_enzyme_ reactions, which are
%       after EX_protein_ but before enzymeForm_
        for k = (length(fullRxn)+length(fullProtein)+1):(length(fullRxn)+length(fullProtein)+length(fullEnzyme))
            if model.C(j,k) ~= 0
                model.C(j,k) = -f * k_eff(k-length(fullRxn)-length(fullProtein));
            end
        end
    end
    
%   Solve and collect results
    FBAsol = optimizeCbModel(model,'max');
    
    result_keff(:,i) = k_eff;
    result_rxn(:,i) = FBAsol.v(1:length(fullRxn));
    result_protein(:,i) = FBAsol.v((length(fullRxn)+1): (length(fullRxn)+length(fullProtein)));
    result_enzyme(:,i) = FBAsol.v((length(fullRxn)+length(fullProtein)+1): (length(fullRxn)+length(fullProtein)+length(fullEnzyme)));
    
%   Progress counter
    if i == noSample*0.1
        fprintf('10%%...');
    elseif i == noSample*0.2
        fprintf('20%%...');
    elseif i == noSample*0.3
        fprintf('30%%...');
    elseif i == noSample*0.4
        fprintf('40%%...');
    elseif i == noSample*0.5
        fprintf('50%%...');
    elseif i == noSample*0.6
        fprintf('60%%...');
    elseif i == noSample*0.7
        fprintf('70%%...');
    elseif i == noSample*0.8
        fprintf('80%%...');
    elseif i == noSample*0.9
        fprintf('90%%...');
    end

end
fprintf('Complete\n');

end