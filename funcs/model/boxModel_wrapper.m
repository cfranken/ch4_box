%<%%                          =======================================================================
%%%                           = boxModel_wrapper.m
%%%			      = Alex Turner
%%%			      = 04/12/2016
%%%			      =----------------------------------------------------------------------
%%%			      = NOTES
%%%			      =  ( 1): A wrapper for the box model.
%%%			      =  ( 2): Precomputes terms for the box model for speed.
%%%			      =  ( 3): Runs the box model.
%%%			      =----------------------------------------------------------------------
%%%			      = INPUTS
%%%			      =  ( 1): St     -- Our time vector.
%%%			      =  ( 2): ems    -- Structure with emission sources (and OH) for the box model.
%%%			      =  ( 3): IC     -- Initial conditions for the box model.
%%%			      =  ( 4): params -- Structure with parameters for the box model.
%%%			      =----------------------------------------------------------------------
%%%			      = OUTPUTS
%%%			      =  ( 1): out -- Structure with the simulated concentrations.
%%%			      =======================================================================

function [ out, S ]	      = boxModel_wrapper(St,S,IC,params)

%%% For the OH feedback
global interactive_OH fixedOH

%%% Set up the emissions for the box model
% Convert CH4, MCF, N2O, C2H6, OH, and CO emissions to units of per day
S(:,[1,2,5,6,7,8,9,10,13,14]) = S(:,[1,2,5,6,7,8,9,10,13,14]) / params.YrToDay;
% 12CH4
S(:,1)			      = 2/params.mm_ch4*S(:,1);    % NH (factor of two is for mass of one hemisphere)
S(:,2)			      = 2/params.mm_ch4*S(:,2);    % SH
% 13CH4
S(:,3)			      = S(:,1).*(1 + S(:,3)/1000); % NH (precompute S1*S2)
S(:,4)			      = S(:,2).*(1 + S(:,4)/1000); % SH (precompute S4*S5)
% MCF
S(:,5)			      = 2/params.mm_mcf*S(:,5);    % NH
S(:,6)			      = 2/params.mm_mcf*S(:,6);    % SH
% N2O
S(:,7)			      = 2/params.mm_n2o*S(:,7);    % NH
S(:,8)			      = 2/params.mm_n2o*S(:,8);    % SH
% C2H6
S(:,9)			      = 2/params.mm_c2h6*S(:,9);   % NH
S(:,10)			      = 2/params.mm_c2h6*S(:,10);  % SH
% OH
if interactive_OH
    S(:,11)		      = 2/params.mm_oh*S(:,11)/params.YrToDay;   % NH
    S(:,12)		      = 2/params.mm_oh*S(:,12)/params.YrToDay;   % SH
else
    S(:,11)		      = params.gmOH*params.DaysToS/params.RxNconv*S(:,11);   % NH
    S(:,12)		      = params.gmOH*params.DaysToS/params.RxNconv*S(:,12);   % SH
end
% CO
S(:,13)			      = 2/params.mm_co*S(:,13);    % NH
S(:,14)			      = 2/params.mm_co*S(:,14);    % SH
% Strat/Trop exchange
tau_TS			      = S(:,15) * params.YrToDay;
% Arbitrary reaction with OH
kX_NH			      = S(:,16) * (60 * 60 * 24);    % NH
kX_SH			      = S(:,17) * (60 * 60 * 24);    % SH



%%% Run the box model with ode45
[T, F]			      = ode15s(@(t,y) boxModel(t,y,St,S,tau_TS,kX_NH,kX_SH,params),params.Tspan,IC,params.odeOpts);

%%% Isotope conversion
F(:,3)			      = (F(:,3) ./F(:,1) -1)*1000;
F(:,4)			      = (F(:,4) ./F(:,2) -1)*1000;
F(:,17)			      = (F(:,17)./F(:,15)-1)*1000;
F(:,18)			      = (F(:,18)./F(:,16)-1)*1000;

%%% Make the output structure
% Do we need to interpolate?
if any(T ~		      = St)
    [pindex,index,slope]      = lininterp1_ind(T,St);
    F			      = F(pindex,:) * (1 - slope) + slope * F(index,:);
end
% Store the result
out.nh_ch4		      = F(:,1);
out.sh_ch4		      = F(:,2);
out.nh_ch4c13		      = F(:,3);
out.sh_ch4c13		      = F(:,4);
out.nh_mcf		      = F(:,5);
out.sh_mcf		      = F(:,6);
out.nh_n2o		      = F(:,7);
out.sh_n2o		      = F(:,8);
out.nh_c2h6		      = F(:,9);
out.sh_c2h6		      = F(:,10);
out.nh_oh		      = F(:,11) * params.n_air/1e9; % Convert to molec/cm3
out.sh_oh		      = F(:,12) * params.n_air/1e9; % Convert to molec/cm3
if ~interactive_OH && ~fixedOH
    out.nh_oh		      = S(:,11) ./ (params.gmOH*params.DaysToS/params.RxNconv);
    out.sh_oh		      = S(:,12) ./ (params.gmOH*params.DaysToS/params.RxNconv);
end
out.nh_co		      = F(:,13);
out.sh_co		      = F(:,14);
out.nh_ch4_strat	      = F(:,15);
out.sh_ch4_strat	      = F(:,16);
out.nh_ch4c13_strat	      = F(:,17);
out.sh_ch4c13_strat	      = F(:,18);
out.nh_mcf_strat	      = F(:,19);
out.sh_mcf_strat	      = F(:,20);
out.nh_n2o_strat	      = F(:,21);
out.sh_n2o_strat	      = F(:,22);
out.nh_c2h6_strat	      = F(:,23);
out.sh_c2h6_strat	      = F(:,24);
out.nh_oh_strat		      = F(:,25) * params.n_air/1e9; % Convert to molec/cm3;
out.sh_oh_strat		      = F(:,26) * params.n_air/1e9; % Convert to molec/cm3;
out.nh_co_strat		      = F(:,27);
out.sh_co_strat		      = F(:,28);

end


%%%			      =======================================================================
%%%			      = END
%%%			      =======================================================================
