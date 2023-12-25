function milp = changeNumKO(milp,numKO,forGRB)

% Update numKO for PC-OptKnock model struct
% 
% USAGE:
% 
%   model_pcok = changeNumKO(milp_pcok,numKO);
%   model_pcok = changeNumKO(milp_pcok,numKO,forGRB);
% 
% INPUTS: 
%   milp:  MILP model struct constructed by formulatePCOptKnock.m
%   numKO: The new maximum number of knocked out allowed
% 
% OPTIONAL INPUTS:
%   forGRB: If the milp model is formulated for Gurobi (LEGACY OPTION)
%           Default: false
% 
% OUTPUTS:
%   milp: MILP model struct with updated numKO
%  
% .. AUTHOR: Herbert Yao, Dec 2023
% 

yIdx = find(startsWith(milp.varnames,'y_'));

if ~exist('forGRB','var')
    milp.b(find(strcmp(milp.constrnames,'max_KO'))) = length(yIdx) - numKO;
elseif ~forGRB
    milp.b(find(strcmp(milp.constrnames,'max_KO'))) = length(yIdx) - numKO;
else
    milp.rhs(find(strcmp(milp.constrnames,'max_KO'))) = length(yIdx) - numKO;
end
end
