%%% =======================================================================
%%% = boxModel.m
%%% = Alex Turner
%%% = 04/12/2016
%%% =----------------------------------------------------------------------
%%% = NOTES
%%% =  ( 1): 2-box model for methane, delta13C, and methylchloroform.
%%% =  ( 2): Adapted from C. Frankenberg's box model.
%%% =----------------------------------------------------------------------
%%% = INPUTS
%%% =  ( 1): t      -- Box model time.
%%% =  ( 2): y      -- Concentrations at time t.
%%% =  ( 3): St     -- Out time vector.
%%% =  ( 4): S      -- Emission sources for the box model.
%%% =  ( 5): tau_ST -- Strat-trop exchange (function of time).
%%% =  ( 6): params -- Structure with the parameters for the box model.
%%% =----------------------------------------------------------------------
%%% = OUTPUTS
%%% =  ( 1): dy -- Changes in the box model concentrations.
%%% =======================================================================

function [ dy ] = boxModel(t,y,St,S,tau_TS,kX_NH,kX_SH,params)

%%% For the OH feedback and the stratosphere
global interactive_OH use_strat fixedOH

%%% Initialize dy for troposphere & stratosphere (strat is same as trop)
% - Column 1:  12CH4 in the Northern Hemisphere
% - Column 2:  12CH4 in the Southern Hemisphere
% - Column 3:  13CH4 in the Northern Hemisphere
% - Column 4:  13CH4 in the Southern Hemisphere
% - Column 5:    MCF in the Northern Hemisphere
% - Column 6:    MCF in the Southern Hemisphere
% - Column 7:    N2O in the Northern Hemisphere
% - Column 8:    N2O in the Southern Hemisphere
% - Column 9:   C2H6 in the Northern Hemisphere
% - Column 10:  C2H6 in the Southern Hemisphere
% - Column 11:    OH in the Northern Hemisphere
% - Column 12:    OH in the Southern Hemisphere
% - Column 13:    CO in the Northern Hemisphere
% - Column 14:    CO in the Southern Hemisphere
dy = zeros(28,1);

%%% Interpolate sources and strat-trop exchange in time
ind = find(t == St);        % Find the index for this timestep
if ~isempty(ind)            % Do we already have it?
    S      = S(ind(1),:);
    tau_TS = tau_TS(ind(1));
    kX_NH  = kX_NH(ind(1));
    kX_SH  = kX_SH(ind(1));
else                        % If not, we'll interpolate
    [pindex,index,slope] = lininterp1_ind(St,t);
    S      = S(pindex,:)      * (1 - slope) + slope * S(index,:);
    tau_TS = tau_TS(pindex,:) * (1 - slope) + slope * tau_TS(index,:);
    kX_NH  = kX_NH(pindex,:) * (1 - slope) + slope * kX_NH(index,:);
    kX_SH  = kX_SH(pindex,:) * (1 - slope) + slope * kX_SH(index,:);
end
% Define the North-South exchange rates
tau_NS       = params.tau_NS;
tau_NS_strat = params.tau_NS_strat;

%%% Conserve mass between strat & trop, so:
% tau_ST = tau_TS/((p_surface-p_tropopause)/(p_tropopause-p_stratopause))
% where m_fac = (p_surface-p_tropopause)/(p_tropopause-p_stratopause)
%             = 5.7047
tau_ST = tau_TS / 5.7047;

%%% Is OH simulated or prescribed?
if ~interactive_OH && ~fixedOH
    y(11) = S(11);
    y(12) = S(12);
end
% Tropospheric reaction rates (y(11) & y(12) are the OH concentrations)
k_12ch4_NH   = (y(11) * params.k_12ch4);   % NH
k_12ch4_SH   = (y(12) * params.k_12ch4);   % SH
k_13ch4_NH   = (y(11) * params.k_13ch4);   % NH
k_13ch4_SH   = (y(12) * params.k_13ch4);   % SH
k_mcf_NH     = (y(11) * params.k_mcf  );   % NH
k_mcf_SH     = (y(12) * params.k_mcf  );   % SH
k_c2h6_NH    = (y(11) * params.k_c2h6 );   % NH
k_c2h6_SH    = (y(12) * params.k_c2h6 );   % SH
k_co_NH      = (y(11) * params.k_co   );   % NH
k_co_SH      = (y(12) * params.k_co   );   % SH
k_ch4_other  = params.k_ch4_other;         % Other methane loss pathways (chlorine in the MBL and soils)
k_c2h6_other = params.k_c2h6_other;        % Other ethane loss pathways
k_co_other   = params.k_co_other;          % Other CO loss pathways
% Stratospheric reaction rates
k_ch4_strat_NH = params.k_ch4_strat_nh;    % NH
k_ch4_strat_SH = params.k_ch4_strat_sh;    % SH
k_mcf_strat    = params.k_mcf_strat;       % Single reaction rate
k_n2o_NH       = params.k_n2o_nh;          % NH
k_n2o_SH       = params.k_n2o_sh;          % SH
k_c2h6_strat   = params.k_c2h6_strat;      % Single reaction rate
k_oh_strat     = params.k_oh_strat;        % Single reaction rate
k_co_strat     = params.k_co_strat;        % Single reaction rate


%%% Compute dy in the troposphere
% 12CH4
dy(1)  = S(1)  + (y(2) - y(1))/tau_NS + (y(15)- y(1))/tau_TS - y(1)*(k_12ch4_NH+k_ch4_other);  % NH
dy(2)  = S(2)  + (y(1) - y(2))/tau_NS + (y(16)- y(2))/tau_TS - y(2)*(k_12ch4_SH+k_ch4_other);  % SH
% 13CH4
dy(3)  = S(3)  + (y(4) - y(3))/tau_NS + (y(17)- y(3))/tau_TS - y(3)*(k_13ch4_NH+k_ch4_other);  % NH
dy(4)  = S(4)  + (y(3) - y(4))/tau_NS + (y(18)- y(4))/tau_TS - y(4)*(k_13ch4_SH+k_ch4_other);  % SH
% MCF
dy(5)  = S(5)  + (y(6) - y(5))/tau_NS + (y(19)- y(5))/tau_TS - y(5)*k_mcf_NH;                  % NH
dy(6)  = S(6)  + (y(5) - y(6))/tau_NS + (y(20)- y(6))/tau_TS - y(6)*k_mcf_SH;                  % SH
% N2O
dy(7)  = S(7)  + (y(8) - y(7))/tau_NS + (y(21)- y(7))/tau_TS;                                  % NH (no loss in troposphere)
dy(8)  = S(8)  + (y(7) - y(8))/tau_NS + (y(22)- y(8))/tau_TS;                                  % SH (no loss in troposphere)
% C2H6
dy(9)  = S(9)  + (y(10)- y(9))/tau_NS + (y(23)- y(9))/tau_TS - y(9)*(k_c2h6_NH+k_c2h6_other);  % NH
dy(10) = S(10) + (y(9) -y(10))/tau_NS + (y(24)-y(10))/tau_TS - y(10)*(k_c2h6_SH+k_c2h6_other); % SH
% OH (allow the feedback?)
if interactive_OH 
dy(11) = S(11)  - y(1)*k_12ch4_NH   - y(13)*k_co_NH - y(11)*kX_NH; % NH
dy(12) = S(12)  - y(2)*k_12ch4_SH   - y(14)*k_co_SH - y(12)*kX_SH; % SH
end
% CO
dy(13) = S(13) + (y(14)-y(13))/tau_NS + (y(27)-y(13))/tau_TS - y(13)*(k_co_NH+k_co_other) + y(1)*k_12ch4_NH + y(3)*k_13ch4_NH;    % NH
dy(14) = S(14) + (y(13)-y(14))/tau_NS + (y(28)-y(14))/tau_TS - y(14)*(k_co_SH+k_co_other) + y(2)*k_12ch4_SH + y(4)*k_13ch4_SH;    % SH

%%% Compute dy in the stratosphere (assuming no fractionation in the stratosphere and that all loss can be representated as first-order)
% Are we using the stratosphere?
if use_strat
% 12CH4
dy(15) = (y(16)-y(15))/tau_NS_strat + (y(1) -y(15))/tau_ST - y(15)*k_ch4_strat_NH;  % NH
dy(16) = (y(15)-y(16))/tau_NS_strat + (y(2) -y(16))/tau_ST - y(16)*k_ch4_strat_SH;  % SH
% 13CH4
dy(17) = (y(18)-y(17))/tau_NS_strat + (y(3) -y(17))/tau_ST - y(17)*k_ch4_strat_NH;  % NH
dy(18) = (y(17)-y(18))/tau_NS_strat + (y(4) -y(18))/tau_ST - y(18)*k_ch4_strat_SH;  % SH
% MCF
dy(19) = (y(20)-y(19))/tau_NS_strat + (y(5) -y(19))/tau_ST - y(19)*k_mcf_strat;     % NH
dy(20) = (y(19)-y(20))/tau_NS_strat + (y(6) -y(20))/tau_ST - y(20)*k_mcf_strat;     % SH
% N2O
dy(21) = (y(22)-y(21))/tau_NS_strat + (y(7) -y(21))/tau_ST - y(21)*k_n2o_NH;        % NH
dy(22) = (y(21)-y(22))/tau_NS_strat + (y(8) -y(22))/tau_ST - y(22)*k_n2o_SH;        % SH
% C2H6
dy(23) = (y(24)-y(23))/tau_NS_strat + (y(9) -y(23))/tau_ST - y(23)*k_c2h6_strat;    % NH
dy(24) = (y(23)-y(24))/tau_NS_strat + (y(10)-y(24))/tau_ST - y(24)*k_c2h6_strat;    % SH
% OH
dy(25) = (y(26)-y(25))/tau_NS_strat + (y(11)-y(25))/tau_ST - y(25)*k_oh_strat;      % NH
dy(26) = (y(25)-y(26))/tau_NS_strat + (y(12)-y(26))/tau_ST - y(26)*k_oh_strat;      % SH
% CO
dy(27) = (y(28)-y(27))/tau_NS_strat + (y(13)-y(27))/tau_ST - y(27)*k_co_strat;      % NH
dy(28) = (y(27)-y(28))/tau_NS_strat + (y(14)-y(28))/tau_ST - y(28)*k_co_strat;      % SH
end

end


%%% =======================================================================
%%% = END
%%% =======================================================================