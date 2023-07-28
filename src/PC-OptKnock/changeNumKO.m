function milp = changeNumKO(milp,numKO)

yIdx = find(startsWith(milp.varnames,'y_'));
milp.rhs(find(strcmp(milp.constrnames,'max_KO'))) = length(yIdx) - numKO;

end