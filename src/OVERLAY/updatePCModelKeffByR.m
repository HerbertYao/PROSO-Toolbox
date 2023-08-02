function model_pc_new = updatePCModelKeffByR(model_pc,rValues)
% This function implements QP solutions to a PC-model.
% It produces a N*1 structure contains context-specific PC-models
% 
% USAGE:
% 
%   model_pc_new = updatePCModelKeffByR(model_pc,rValues);
% 
% INPUTS:
% 
%   model_pc: A PC model built from the function 'pcModel.m'
%   rValues:  The r vector that is potentially produced by nonconvex QP
% 
% OUTPUTS:
% 
%   model_pc_new: PC-model with updated R
% 

model_pc_new = model_pc;

fullCplxIdx = find(startsWith(model_pc.mets,'cplx_'));
enzymeFormIdx = find(startsWith(model_pc.rxns,'enzymeForm_'));

if length(fullCplxIdx) ~= length(rValues)
    error('Dimension of r must equal to the number of cplx in the model');
end

for i = 1:length(rValues)
    model_pc_new.S(fullCplxIdx(i),enzymeFormIdx)...
        = model_pc.S(fullCplxIdx(i),enzymeFormIdx) / rValues(i);
end

end