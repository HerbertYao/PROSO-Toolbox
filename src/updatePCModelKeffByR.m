function model_pc_new = updatePCModelKeffByR(model_pc,rValues)

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