function milp = changeNumKO(milp,numKO)
% Update numKO for PC-OptKnock model struct
% 
% USAGE:
% 
%   milp_pcok = changeNumKO(milp_pcok,2);
% 
% INPUTS:
% 
%   model:  MILP model struct constructed by formulatePCOptKnock.m
%   numKO:  The new maximum number of knocked out allowed
% 
% OUTPUTS:
% 
%   milp: MILP model struct with updated numKO
% 

yIdx = find(startsWith(milp.varnames,'y_'));
milp.rhs(find(strcmp(milp.constrnames,'max_KO'))) = length(yIdx) - numKO;

end