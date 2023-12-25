function plotPCFVASolutions(FVAsols,model,rxnList,PltDim,FVAsols_enz,clr)

% Plot PC-FVA results in a compact and easy-to-understand form
% 
% USAGE:
% 
%   plotPCFVASolutions(FVAsols,model_ori,{'EX_glc__D_e','EX_o2_e'},...
%       [1,2],FVAsols_enz)
% 
% INPUTS:
%   FVAsols: 3d matrix as an output from contextSpecificPCFVA.m
%   model:   Unmodified M-model struct or PC-model struct
%   rxnList: List of metabolic rxns to be plotted
%   PltDim:  The dimension of subplots
% 
% OPTIONAL INPUTS:
%   FVAsols_enz: FVA solution of reaction-specific enzyme level. This will
%                be plotted as black lines. Please consult the tutorial
%                about its format. Default = []
%   clr:         User specified color of the bar plot. Must be in a N*1
%                cell with each entry containing a color. It has a default.
% 
% OUTPUTS:
%   N/A
% 
% .. AUTHOR: Herbert Yao, Dec 2023
% 

if ~exist('clr','var')
    clr = {[0.4660 0.6740 0.1880],...
        [0.9290 0.6940 0.1250],...
        [0.8500 0.3250 0.0980],...
        [1 0 0]};
end

recWids = [0.8,0.65,0.5,0.35];

[~,noExp,noOptPerc] = size(FVAsols);
noOptPerc = noOptPerc/2;

figure;

for i = 1:length(rxnList)
    subplot(PltDim(1),PltDim(2),i);
    yyaxis left;
    xlabel('Exp No.');
    ylabel('flux');

    idx = find(strcmp(model.rxns,rxnList{i}));
    
    hold on;
    for j = 1:noExp
        for k = 1:noOptPerc
            rectangle('Position',...
                [(j-recWids(k)/2),...
                FVAsols(idx,j,k),...
                recWids(k),...
                (FVAsols(idx,j,2*noOptPerc+1-k)-FVAsols(idx,j,k))],...
                'FaceColor',clr{k});
        end
    end
    hold off;
    ylabel('flux');

    if ~isempty(model.rules{idx}) && exist('FVAsols_enz','var')
        yyaxis right;
        ylabel('max enzyme abund.');
        p = plot(FVAsols_enz(idx,:),'-k');
        p.LineWidth = 2;
        ylabel('expression');
    end

    xticks(1:noExp);
    title(rxnList{i},'FontSize',10,'Interpreter','none');
end

end