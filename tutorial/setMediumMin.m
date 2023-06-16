function model = setMediumMin(model, identifier)

% Return a COBRA model with bounds properly set to Minimal medium without
% carbon substrates. The resulting growth rate should be zero.
% 
% USAGE:
%   model = setMediumMin(model, 'b')
% 
% INPUTS:
%   model:  COBRA model 1x1 structure
%   identifier: the chemical and reaction id used in the model
%               's' for Seed id
%               'b' for Bigg id
% 
% OUTPUTS:
%   model:  COBRA model set to minimal medium
% 

if strcmp(identifier, 's')
    model = changeRxnBounds(model, 'EX_cpd00001(e)', -1000, 'l');
    model = changeRxnBounds(model, 'EX_cpd00009(e)', -1000, 'l');
    model = changeRxnBounds(model, 'EX_cpd00011(e)', -1000, 'l');
    model = changeRxnBounds(model, 'EX_cpd00021(e)', -1000, 'l');
    model = changeRxnBounds(model, 'EX_cpd00030(e)', -1000, 'l');
    model = changeRxnBounds(model, 'EX_cpd00034(e)', -1000, 'l');
    model = changeRxnBounds(model, 'EX_cpd00048(e)', -1000, 'l');
    model = changeRxnBounds(model, 'EX_cpd00058(e)', -1000, 'l');
    model = changeRxnBounds(model, 'EX_cpd00067(e)', -1000, 'l');
    model = changeRxnBounds(model, 'EX_cpd00149(e)', -1000, 'l');
    model = changeRxnBounds(model, 'EX_cpd00205(e)', -1000, 'l');
    model = changeRxnBounds(model, 'EX_cpd00254(e)', -1000, 'l');
    model = changeRxnBounds(model, 'EX_cpd00528(e)', -1000, 'l');
    model = changeRxnBounds(model, 'EX_cpd00971(e)', -1000, 'l');
    model = changeRxnBounds(model, 'EX_cpd00013(e)', -1000, 'l');
    model = changeRxnBounds(model, 'EX_cpd01012(e)', -1000, 'l');
    model = changeRxnBounds(model, 'EX_cpd10516(e)', -1000, 'l');
    model = changeRxnBounds(model, 'EX_cpd00244(e)', -1000, 'l');
    model = changeRxnBounds(model, 'EX_cpd00007(e)', -20, 'l');

    model = changeRxnBounds(model, 'EX_cpd00001(e)', 1000, 'u');
    model = changeRxnBounds(model, 'EX_cpd00009(e)', 1000, 'u');
    model = changeRxnBounds(model, 'EX_cpd00011(e)', 1000, 'u');
    model = changeRxnBounds(model, 'EX_cpd00021(e)', 1000, 'u');
    model = changeRxnBounds(model, 'EX_cpd00030(e)', 1000, 'u');
    model = changeRxnBounds(model, 'EX_cpd00034(e)', 1000, 'u');
    model = changeRxnBounds(model, 'EX_cpd00048(e)', 1000, 'u');
    model = changeRxnBounds(model, 'EX_cpd00058(e)', 1000, 'u');
    model = changeRxnBounds(model, 'EX_cpd00067(e)', 1000, 'u');
    model = changeRxnBounds(model, 'EX_cpd00149(e)', 1000, 'u');
    model = changeRxnBounds(model, 'EX_cpd00205(e)', 1000, 'u');
    model = changeRxnBounds(model, 'EX_cpd00254(e)', 1000, 'u');
    model = changeRxnBounds(model, 'EX_cpd00528(e)', 1000, 'u');
    model = changeRxnBounds(model, 'EX_cpd00971(e)', 1000, 'u');
    model = changeRxnBounds(model, 'EX_cpd00013(e)', 1000, 'u');
    model = changeRxnBounds(model, 'EX_cpd01012(e)', 1000, 'u');
    model = changeRxnBounds(model, 'EX_cpd10516(e)', 1000, 'u');
    model = changeRxnBounds(model, 'EX_cpd00244(e)', 1000, 'u');
    model = changeRxnBounds(model, 'EX_cpd00007(e)', 1000, 'u');

    model = changeRxnBounds(model, 'EX_cpd00027(e)', 0, 'l');
    model = changeRxnBounds(model, 'EX_cpd00221(e)', 0, 'l');
    model = changeRxnBounds(model, 'EX_cpd00029(e)', 0, 'l');
    model = changeRxnBounds(model, 'EX_cpd00100(e)', 0, 'l');
    model = changeRxnBounds(model, 'EX_cpd00033(e)', 0, 'l');
    model = changeRxnBounds(model, 'EX_cpd00035(e)', 0, 'l');
    model = changeRxnBounds(model, 'EX_cpd00051(e)', 0, 'l');
    model = changeRxnBounds(model, 'EX_cpd00041(e)', 0, 'l');
    model = changeRxnBounds(model, 'EX_cpd00023(e)', 0, 'l');
    model = changeRxnBounds(model, 'EX_cpd00119(e)', 0, 'l');
    model = changeRxnBounds(model, 'EX_cpd00322(e)', 0, 'l');
    model = changeRxnBounds(model, 'EX_cpd00107(e)', 0, 'l');
    model = changeRxnBounds(model, 'EX_cpd00027(e)', 0, 'l');
    model = changeRxnBounds(model, 'EX_cpd00156(e)', 0, 'l');
    model = changeRxnBounds(model, 'EX_cpd00064(e)', 0, 'l');
    model = changeRxnBounds(model, 'EX_cpd00039(e)', 0, 'l');
    model = changeRxnBounds(model, 'EX_cpd00066(e)', 0, 'l');
    model = changeRxnBounds(model, 'EX_cpd00129(e)', 0, 'l');
    model = changeRxnBounds(model, 'EX_cpd00054(e)', 0, 'l');
    model = changeRxnBounds(model, 'EX_cpd00161(e)', 0, 'l');
    model = changeRxnBounds(model, 'EX_cpd00069(e)', 0, 'l');
    
elseif strcmp(identifier, 'b')
    model = changeRxnBounds(model, 'EX_h2o_e', -1000, 'l');
    model = changeRxnBounds(model, 'EX_pi_e', -1000, 'l');
    model = changeRxnBounds(model, 'EX_co2_e', -1000, 'l');
    model = changeRxnBounds(model, 'EX_fe2_e', -1000, 'l');
    model = changeRxnBounds(model, 'EX_mn2_e', -1000, 'l');
    model = changeRxnBounds(model, 'EX_zn2_e', -1000, 'l');
    model = changeRxnBounds(model, 'EX_so4_e', -1000, 'l');
    model = changeRxnBounds(model, 'EX_cu2_e', -1000, 'l');
    model = changeRxnBounds(model, 'EX_h_e', -1000, 'l');
    model = changeRxnBounds(model, 'EX_cobalt2_e', -1000, 'l');
    model = changeRxnBounds(model, 'EX_k_e', -1000, 'l');
    model = changeRxnBounds(model, 'EX_mg2_e', -1000, 'l');
    %model = changeRxnBounds(model, 'EX_n2_e', -1000, 'l'); % Note there are no n2 releasing rxn in Sanjeev's model
    model = changeRxnBounds(model, 'EX_na1_e', -1000, 'l');
    model = changeRxnBounds(model, 'EX_nh4_e', -1000, 'l');
    model = changeRxnBounds(model, 'EX_cd2_e', -1000, 'l');
    model = changeRxnBounds(model, 'EX_fe3_e', -1000, 'l');
    model = changeRxnBounds(model, 'EX_ni2_e', -1000, 'l');
    model = changeRxnBounds(model, 'EX_o2_e', -20, 'l');
    
    model = changeRxnBounds(model, 'EX_glc__D_e', 0, 'l');
    model = changeRxnBounds(model, 'EX_lac__D_e', 0, 'l');
    model = changeRxnBounds(model, 'EX_ac_e', 0, 'l');
    model = changeRxnBounds(model, 'EX_glyc_e', 0, 'l');
    model = changeRxnBounds(model, 'EX_gly_e', 0, 'l');
    model = changeRxnBounds(model, 'EX_ala__L_e', 0, 'l');
    model = changeRxnBounds(model, 'EX_arg__L_e', 0, 'l');
    model = changeRxnBounds(model, 'EX_asp__L_e', 0, 'l');
    model = changeRxnBounds(model, 'EX_glu__L_e', 0, 'l');
    model = changeRxnBounds(model, 'EX_his__L_e', 0, 'l');
    model = changeRxnBounds(model, 'EX_ile__L_e', 0, 'l');
    model = changeRxnBounds(model, 'EX_leu__L_e', 0, 'l');
    model = changeRxnBounds(model, 'EX_glc__D_e', 0, 'l');
    model = changeRxnBounds(model, 'EX_val__L_e', 0, 'l');
    model = changeRxnBounds(model, 'EX_orn_e', 0, 'l');
    model = changeRxnBounds(model, 'EX_lys__L_e', 0, 'l');
    model = changeRxnBounds(model, 'EX_phe__L_e', 0, 'l');
    model = changeRxnBounds(model, 'EX_pro__L_e', 0, 'l');
    model = changeRxnBounds(model, 'EX_ser__L_e', 0, 'l');
    model = changeRxnBounds(model, 'EX_thr__L_e', 0, 'l');
    model = changeRxnBounds(model, 'EX_tyr__L_e', 0, 'l');
    model = changeRxnBounds(model, 'EX_csn_e', 0, 'l');
    model = changeRxnBounds(model, 'EX_ura_e', 0, 'l');
    
else
    error('Unrecognized identifier');
end

end
