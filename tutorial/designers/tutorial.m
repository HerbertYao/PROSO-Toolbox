% Tutorial for Cell Designer Functions

% This script will provide a detailed guide on how to use functions in 
% OVERLAY/src/Strain-Designer as well as the concepts behind algorithms

initCobraToolbox(false);
changeCobraSolver('gurobi','all',1);
clc;

%% Prepare PC-model

% In this tutorial, we will use S. cerevisiae metabolic reconstruction 
% iMM904 as an example 

if ~exist('model_ori','var')
%     model_ori = readCbModel('iML1515.xml');
    model_ori = readCbModel('iMM904.xml');
end

% rxnIdx.biomass = 2669;
% rxnIdx.EX_glc = 181;
% rxnIdx.EX_o2 = 1982;
% rxnIdx.EX_succ = 594;

rxnIdx.biomass = 1521;
rxnIdx.EX_glc = 508;
rxnIdx.EX_o2 = 599;
rxnIdx.EX_succ = 594;

% This model is setup to solely grow on glucose

% model_ori.lb(rxnIdx.EX_o2) = -1000; % alleviate oxygen constraint

FBAsol = optimizeCbModel(model_ori,'max');
fprintf('With glucose: growth rate: %.3f hr^-1\n',FBAsol.f);
fprintf('              glucose consumption rate: %.3f mmol/gDW/h\n',FBAsol.v(rxnIdx.EX_glc));
fprintf('              oxygen consumption rate: %.3f mmol/gDW/h\n',FBAsol.v(rxnIdx.EX_o2));

% Check if glucose is truly the sole carbon source

model = model_ori;
model.lb(rxnIdx.EX_glc) = 0;

FBAsol = optimizeCbModel(model,'max');
if ~FBAsol.stat
    fprintf('Without glucose: infeasible\n');
else
    fprintf('Without glucose: growth rate: %.3f hr^-1\n',FBAsol.f);
    fprintf('                 glucose consumption rate: %.3f mmol/gDW/h\n',FBAsol.v(rxnIdx.EX_glc));
    fprintf('                 oxygen consumption rate: %.3f mmol/gDW/h\n',FBAsol.v(rxnIdx.EX_o2));
end

% This is not an demonstration of OVERLAY, so PC-model formulation will be
% greatly simplified. 
% 
% The 'pcModel' function can also formulate PC-model without specified 
% protein FASTA file, but instead provide an average protein length for 
% protein molecular weight calculation. This will reduce the rigor of the 
% resulted PC-model, but it can be used for testing the workflow and
% investigating the feasibility of study before starting. Here I will
% utilize this feature.

if ~exist('model_pc_ori','var')
    model_pc_ori = pcModel(model_ori,200,150);
end

% Here we specify the average length of proteins in E. coli to be 200 aa
% long, which substitutes for the FASTA filename. The protein weight
% constraint is set at 150mg/gDW or 15% of dry weight. Details about this 
% function are explained in the other tutorial.

% We will not refine the PC-model here because this is only a showcase for
% cell designer functions. 

%% PC-OptKnock

% PC-OptKnock has identical idea to OptKnock but deals with PC-model. This
% allows more practical gene-knockout suggestions from the solver. 

% Functions in this section can be very time consuming to run and solve, so
% please have one of your better device and take your time. 

% The model struct for OptKnock or PC-OptKnock can be very handily 
% constructed by calling a single function. Let's try doing the regular
% OptKnock first by using M-model and specifying the PCBool 
% (PCBool = false). 

if ~exist('model_ok','var')
    [model_ok,param] = formulatePCOptKnock(model_ori,1,0.05,true,false);
end

% Let's inspect the MILP model's structure:
% 
% Variables (15437 in total):
%   - 2712 metabolic reactions
%   - 1877 w variables, one for each metabolite
%   - 2712 a, b, y, and sigma variables, respectively, one set of these 
%     variables for each metabolic reaction
% 
% Constraints (20863 in total):
%   - 1877 metabolites
%   - 2712 dual constraints, one for each metabolic reaction
%   - 2712 pairs of eta and mu constraints, one pair for each metabolic
%     reactions
%   - 2712 pairs of sigma_L and sigma_U constraints, one pair for each
%     metabolic reactions
%   - 1 strong duality constraint
%   - 2712 pairs of lb_y and ub_y constraints, one pair for each metabolic
%     reactions
%   - 1 maxKO constraint

% Before solving this MILP, we still need to set the objective function

if ~exist('OKsol','var')
    model_ok.obj(rxnIdx.EX_succ) = 1;
    model_ok = changeNumKO(model_ok,1);
    OKsol = gurobi(model_ok,param);
end

% Now find knocked out reactions and implement it

idx = find(~OKsol.x(find(startsWith(model_ok.varnames,'y_'))));
model_test = model_ori;
model_test.lb(idx) = 0;
model_test.ub(idx) = 0;

drawSolutionEnvelope(model_test,model_ori.rxns{rxnIdx.biomass},model_ori.rxns{rxnIdx.EX_succ},20,true);

% This is equivalent to the original OptKnock

% =========================================================================
% 
% Now let's try PC-OptKnock
% 
% Construct the MILP struct for PC-OptKnock from PC-model while setting the
% PCBool to true. This will take significantly longer than constructing the
% MILP struct for the regular OptKnock. 

if ~exist('model_pcok','var')
    [model_pcok,param] = formulatePCOptKnock(model_pc_ori,1,0.05,true,true);
end

yIdx = find(startsWith(model_pcok.varnames,'y_')); % binary var index
models_sol = {};

% Let's now first try doing only 1 knockout. You ALWAYS want to start with
% K = 1 and work your way up. 

K_init = 1;

model_pcok.obj(rxnIdx.EX_succ) = 1; % specify the outer objective
model_pcok = changeNumKO(model_pcok,K_init); % set K
PCOKsol = gurobi(model_pcok,param);

% If we do get a growth-coupling solution (we should according to the
% design), we can get a copy of PC-model and implement the knockout to see
% if it is indeed growth-coupled. 

if PCOKsol.x(rxnIdx.EX_succ) > 1e-6

    idx = find(~PCOKsol.x(yIdx));
    models_sol{1} = model_pc_ori;
    models_sol{1}.lb(idx) = 0;
    models_sol{1}.ub(idx) = 0;
    
    FBAsol = optimizeCbModel(models_sol{1},'max');
    fprintf('K = %d: EX_succ_e = %.2f, mu = %.2f\n',K_init,FBAsol.v(rxnIdx.EX_succ),FBAsol.v(rxnIdx.biomass));
    fprintf('Closed proteins: %s\n',model_pc_ori.rxns{idx(1)});
    drawSolutionEnvelope(models_sol{1},model_pc_ori.rxns{rxnIdx.biomass},model_pc_ori.rxns{rxnIdx.EX_succ},20,true);
else
    error('No growth coupling solution found by PC-OptKnock. Please consider increase K');
end

% In this case, we find a growth-coupling solution with only K = 1, which 
% is an ideal scenario. If PCOKsol.x(rxnIdx.EX_succ) is almost 0, you should
% increase K by 1 and re-run this section. 
% 
% Since we have found a solution, it is attempting to increase the
% productivity. The most computationally efficient way is to set K = 2 and
% fixing the previous growth-coupling solution. We call this trial K = 1+1.
% This finds a local optima starting from the existing solution of K = 1, 
% which is very likely less optimized than computing the true K = 2 but 
% much more efficient. 
% 
% So let's do K = 1+1 here

model_pcok_alt = changeNumKO(model_pcok,K_init+1);
model_pcok_alt.ub(yIdx(idx)) = 0; % forcing the previous knockout
PCOKsol = gurobi(model_pcok_alt,param);

idx = find(~PCOKsol.x(find(startsWith(model_pcok_alt.varnames,'y_'))));
models_sol{2} = model_pc_ori;
models_sol{2}.lb(idx) = 0;
models_sol{2}.ub(idx) = 0;

FBAsol = optimizeCbModel(models_sol{2},'max');
fprintf('K = %d+1: EX_succ_e = %.2f, mu = %.2f\n',K_init,FBAsol.v(rxnIdx.EX_succ),FBAsol.v(rxnIdx.biomass));
fprintf('Closed proteins: %s, %s\n',model_pc_ori.rxns{idx(1)},model_pc_ori.rxns{idx(2)});
drawSolutionEnvelope(models_sol{2},model_pc_ori.rxns{rxnIdx.biomass},model_pc_ori.rxns{rxnIdx.EX_succ},20,true);

% And then K = 1+1+1: 

model_pcok_alt = changeNumKO(model_pcok,K_init+2);
model_pcok_alt.ub(yIdx(idx)) = 0;
PCOKsol = gurobi(model_pcok_alt,param);

idx = find(~PCOKsol.x(find(startsWith(model_pcok_alt.varnames,'y_'))));
models_sol{3} = model_pc_ori;
models_sol{3}.lb(idx) = 0;
models_sol{3}.ub(idx) = 0;

FBAsol = optimizeCbModel(models_sol{3},'max');
fprintf('K = %d+1+1: EX_succ_e = %.2f, mu = %.2f\n',K_init,FBAsol.v(rxnIdx.EX_succ),FBAsol.v(rxnIdx.biomass));
fprintf('Closed proteins: %s, %s, %s\n',model_pc_ori.rxns{idx(1)},model_pc_ori.rxns{idx(2)},model_pc_ori.rxns{idx(3)});
drawSolutionEnvelope(models_sol{3},model_pc_ori.rxns{rxnIdx.biomass},model_pc_ori.rxns{rxnIdx.EX_succ},20,true);

% And then K = 1+1+1+1: 

model_pcok_alt = changeNumKO(model_pcok,K_init+3);
model_pcok_alt.ub(yIdx(idx)) = 0;
PCOKsol = gurobi(model_pcok_alt,param);

idx = find(~PCOKsol.x(find(startsWith(model_pcok_alt.varnames,'y_'))));
models_sol{4} = model_pc_ori;
models_sol{4}.lb(idx) = 0;
models_sol{4}.ub(idx) = 0;

FBAsol = optimizeCbModel(models_sol{4},'max');
fprintf('K = %d+1+1+1: EX_succ_e = %.2f, mu = %.2f\n',K_init,FBAsol.v(rxnIdx.EX_succ),FBAsol.v(rxnIdx.biomass));
fprintf('Closed proteins: %s, %s, %s, %s\n',...
    model_pc_ori.rxns{idx(1)},model_pc_ori.rxns{idx(2)},model_pc_ori.rxns{idx(3)},model_pc_ori.rxns{idx(4)});
drawSolutionEnvelope(models_sol{4},model_pc_ori.rxns{rxnIdx.biomass},model_pc_ori.rxns{rxnIdx.EX_succ},20,true);

% As you can see, we could keep stacking up with minor computational effort
% and achieve much better growth-coupling strain. 
% If we try to solve K = 4 from the beginning, it will cost more than 100x
% resources, although it can guarantee to reach a global optima, of which
% in this case we don't neccesarily need to achieve anyway.

%% MOPA

% The algorithm of MOPA (minimization of proteome adjustment) is adapted 
% from MOMA (minimization of metabolic adjustment). This suggests that once
% a genomic modification is made, the proteome (metabolism) will gradually
% shift from the old operating point toward the new. As a consequence, the
% short-term behavior of the modified strain may appear to be sub-optimal,
% because the cell is still operating around the old operating point. 
% 
% Let's use the PC-OptKnock K = 1+1+1 result as an example, note that these
% are the knockouts from that tutorial section:

KOProteins = {'EX_protein_YPR002W','EX_protein_Q0250','EX_protein_YBR196C','EX_protein_YKL029C'};
% KOProteins = {'EX_protein_Q0250','EX_protein_YBR196C'};

proteinExIdx = find(startsWith(model_pc_ori.rxns,'EX_protein_'));
metRxnIdx = 1:1577;

% MOPA gives the MOPA solution, wildtype optimum, and KO strain optimum

[MOPAsol,FBAsol_wt,FBAsol_mut] = MOPA(model_pc_ori,KOProteins);

% It is difficult to understand what's going on simply by looking at this
% solution vectors. Here we can use a tweening algorithm, which gives the
% optimal route for the cell to progressively move its proteome toward the
% new operating point. Note in the proteinTween algorithm, the protein is
% re-allocated uniformly between each step, but the metabolism change 
% usually appear non-uniformal. 

FBAsols = proteinTween(models_sol{4},MOPAsol,'BIOMASS_SC5_notrace',20);

% The solution is too high in dimension to visualize, but we can plot the 
% manifold after PCA. Note that PCA need to be done for metabolism and
% proteome separately.

[coef_m,sc_m,~,~,exp_m] = pca([FBAsol_wt.v(metRxnIdx),FBAsols(metRxnIdx,:)]');
[coef_p,sc_p,~,~,exp_p] = pca([FBAsol_wt.v(proteinExIdx),FBAsols(proteinExIdx,:)]');

% Plot the manifold

figure;

subplot(2,1,1);
hold on;
plot(sc_m(2:end,1),sc_m(2:end,2),'o-');
plot(sc_m(1,1),sc_m(1,2),'*','MarkerSize',10);
plot(sc_m(end,1),sc_m(end,2),'*','MarkerSize',10);
plot(sc_m(2,1),sc_m(2,2),'*','MarkerSize',10);
legend({'migrating route','wildtype optimum point','knockout strain optimum point','MOPA point'});
title('Metabolism');
xlabel(['PC 1: ',num2str(exp_m(1)),'%']);
ylabel(['PC 2: ',num2str(exp_m(2)),'%']);

subplot(2,1,2);
hold on;
plot(sc_p(2:end,1),sc_p(2:end,2),'o-');
plot(sc_p(1,1),sc_p(1,2),'*','MarkerSize',10);
plot(sc_p(end,1),sc_p(end,2),'*','MarkerSize',10);
plot(sc_p(2,1),sc_p(2,2),'*','MarkerSize',10);
legend({'migrating route','wildtype optimum point','knockout strain optimum point','MOPA point'});
title('Proteome');
xlabel(['PC 1: ',num2str(exp_p(1)),'%']);
ylabel(['PC 2: ',num2str(exp_p(2)),'%']);

% From this plot, you can indeed see the MOPA solution is the closest to
% the wildtype point, but only for the proteome. The pattern is a lot
% harder to see for metabolic vectors. However, the re-allocation of
% proteome, despite not moving linearly itself, allows an almost linear 
% movement for the metabolism toward the new optimum. This is what we want
% to achieve by MOMA.

% We can also see how growth and succinate production are affected

figure;

hold on;
plot3(FBAsols(rxnIdx.biomass,:),sc_p(2:end,1),FBAsols(rxnIdx.EX_succ,:),'o-');
plot3(FBAsol_wt.v(rxnIdx.biomass),sc_p(1,1),FBAsol_wt.v(rxnIdx.EX_succ),'*','MarkerSize',10);
plot3(FBAsols(rxnIdx.biomass,end),sc_p(end,1),FBAsols(rxnIdx.EX_succ,end),'*','MarkerSize',10);
plot3(FBAsols(rxnIdx.biomass,2),sc_p(2,1),FBAsols(rxnIdx.EX_succ,2),'*','MarkerSize',10);
hold off;
legend({'migrating route','wildtype optimum point','knockout strain optimum point','MOPA point'});
xlabel('Growth Rate');
ylabel('Proteome PC 1');
zlabel('Succinate Production Rate');
% set(gca,'Xscale','log');

%% minimalGenome

% The minimal genome is an interesting concept in synthetic biology: by
% reducing the size of the genome while preserving the desired synthetic
% pathway, we can largely avoid regulatory complications and achieve and
% maintain a relatively high productivity. 
% 
% The in-silico minimal genome simulation can be easily done as a MILP
% problem from a PC-model. We are aiming to shutdown as many proteins while
% keeping the growth rate above a threshold, which is conceptually simple.
% However, the MILP can become very expensive to solve as we trying to
% knockout too many proteins at once. Thus, we can also take the step-wise
% approach, although this does not gaurantee to find the global
% optimality, just like the step-wise PC-OptKnock didn't. 
% 
% Let's use the wildtype S. cerevisiae as an example. We want to reduce its
% genome size while preserving the ability to produce succinate and grow to
% a certain degree. 

