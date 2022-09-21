function [model_ok,koList] = pcOptKnock(pcOptKnockMILP,objIdx,bmIdx,bmMin,K,protectList)

% Finding a growth coupled mutant strain of a given milp by knocking
% out a designated number of proteins. Input milp must be prepared by the
% function formulatePCOptKnock.m
% 
% USAGE:
% 
%   model_min = findMinimalCell(model_pc,0.5);
% 
% INPUTS:
% 
%   model:       A PC-model produced by function pcModel.m or preferably
%                refined by adjustStoichAndKeff.m
%   objRxn:      
%   bmRxn:
%   K:           The number of protein knock-out allowed so must be an
%                positive integer
%   protectList: List of proteins being protected from knocking out
% 
% OUTPUTS:
% 
%   model_ok:
%   KOList:
% 
% Note: This function can be very computationally expensive especially with
% large K values on a large pc-model. Please try to start from K=1 and
% working the way up if a solution isn't found.
% 



end