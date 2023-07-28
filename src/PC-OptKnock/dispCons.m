function dispCons(model,n)

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