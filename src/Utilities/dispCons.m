function dispCons(model,n)
% A utility function for inspecting the n-th constraint of a MILP model.
% This meant to be a MILP version of surfNet() function in COBRA Toolbox
% 
% USAGE:
% 
%   dispCons(model_pcok,8000);
% 
% INPUTS:
% 
%   model: A MILP model formulated by other functions in this toolbox
%   n:     The index of constraint to be displayed.
% 
% OUTPUTS:
% 
%   Constraint dual_cplxForm_x(514): 
%     -1*w_protein_YLL062C + w_cplx_x(514) + -1*a_cplxForm_x(514) 
%     + b_cplxForm_x(514) + sigma_cplxForm_x(514) 
%     = 0
% 

fprintf('Constraint %s: \n  ',model.constrnames{n});

idx = find(model.A(n,:));

for i = 1:length(idx)
    if model.A(n,idx(i)) ~= 1
        fprintf('%d*',full(model.A(n,idx(i))));
    end
    fprintf('%s ',model.varnames{idx(i)});

    if i ~= length(idx)
        fprintf('+ ');
    end
end

fprintf('\n  ');
fprintf(model.sense(n));
fprintf(' %d\n',model.rhs(n))

end