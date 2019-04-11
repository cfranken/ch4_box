%%% =======================================================================
%%% = getParameters.m
%%% = Alex Turner
%%% = 04/12/2016
%%% =----------------------------------------------------------------------
%%% = NOTES
%%% =  ( 1): Computes the parameters that are needed for the box model.
%%% =----------------------------------------------------------------------
%%% = INPUTS
%%% =  ( 1): St -- Our time vector.
%%% =----------------------------------------------------------------------
%%% = OUTPUTS
%%% =  ( 1): params -- Structure with the parameters for the box model.
%%% =======================================================================

function [ params ] = getParameters(St)

%%% Are we using non-OH sinks?
global use_other_sinks

%%% Box model parameters and unit conversions
% Unit conversions
DaysToS = 60 * 60 * 24;         % Days to Seconds
YrToDay = 365.25;               % Years to Days
% Masses
m       = 5.15e21;              % Rough: Total mass of atmosphere in g
n_air   = 2.69e19;              % Rough: number density of air (molec/cm3)
m_air   = 28.8;                 % Average molar mass of the atmosphere in g/mol (mostly N2 and O2)
m_ch4   = 16;                   % Molar mass of CH4
m_mcf   = 133.4;                % Molar mass of methyl-chloroform
m_n2o   = 44.013;               % Molar mass of nitrous oxide
m_c2h6  = 30.07;                % Molar mass of ethane
m_oh    = 17.008;               % Molar mass of hydroxyl
m_co    = 28.01;                % Molar mass of carbon monoxide
mConv   = m/m_air/1e12*1e-9;	% Factor to convert mixing ratios to Tg
mm_ch4  = m_ch4  * mConv;       % Convert CH4 mixing ratios (ppb) to Tg
mm_mcf  = m_mcf  * mConv;       % Convert MCF mixing ratios (ppt) to Gg
mm_n2o  = m_n2o  * mConv;       % Convert N2O mixing ratios (ppb) to Tg
mm_c2h6 = m_c2h6 * mConv;       % Convert C2H6 mixing ratios (ppt) to Gg
mm_oh   = m_oh  * mConv;        % Convert OH number densities (ppb) to Tg
mm_co   = m_co  * mConv;        % Convert CO mixing ratios (ppb) to Tg
% Kinetic isotope effect
KIE  = 1.005;                   % Kinetic isotope effect (k12/k13) from Burkholder et al. (2015)
iKIE = 1/KIE;                   % Inverse of KIE
% Gloal mean OH
gmOH = 1.0d6;                   % Global mean OH (molec/cm^3)
% Reaction rate conversion factor
RxNconv = n_air / 1d9;          % Convert ppb to molec/cm3
% Reaction rates and lifetimes
k_12ch4        = 3.395e-15;     % reaction rate of OH with 12CH4 at 270K (cm3/molec/s)
k_13ch4        = k_12ch4 * iKIE;% reaction rate of OH with 13CH4 (cm3/molec/s)
t_n2oNH        = 163;           % N2O lifetime (yr) in the NH (Brown et al., 2013; Table 4)
t_n2oSH        = 109;           % N2O lifetime (yr) in the SH (Brown et al., 2013; Table 4)
t_ch4_strat_nh = 188;           % CH4 lifetime (yr) in the NH (Brown et al., 2013; Table 4)
t_ch4_strat_sh = 200;           % CH4 lifetime (yr) in the SH (Brown et al., 2013; Table 4)
t_mcf_strat    = 31;            % 15% of loss is in stratosphere, using MCF lifetime of 5.5 yr (ratio of losses is proportional to inverse ratio of lifetimes Makide & Rowland, 1981)
k_c2h6         = 1.752e-13;     % reaction rate of OH with C2H6 at 270K (cm3/molec/s)
k_co_A         = 1.0133e-12;    % reaction rate of OH with CO at 270K (cm3/molec/s), sum of bi- & termolecular reactions (Burkholder, Table 2.1)
k_co_B         = 2e-13;         % reaction rate (cm3/molec/s) from Prather (1993)
k_c2h6_other   = 0.15/0.25;     % (yr^-1) Other ethane losses (85% relative to OH loss)
k_c2h6_strat   = 12;            % (yr^-1) Assuming lifetime of 1 month in the stratosphere
k_co_strat     = 12;            % (yr^-1) Assuming lifetime of 1 month in the stratosphere
k_oh_strat     = 12;            % (yr^-1) Assuming lifetime of 1 month in the stratosphere
k_ch4_other    = 1/151;         % (yr^-1) Soil uptake and tropospheric chlorine (6% of total loss; Kirschke et al., 2013)
k_co_other     = 0;             % (yr^-1) Other CO losses (currently neglecting)
k_mcf_A        = 5.66e-15;      % reaction rate with MCF (computed such that lifetime is about 5.5 years, Talukdar et al 1992, taken at 273K)
k_mcf_B        = 6.05e-15;      % reaction rate determined by AJT
% Lifetimes
tau_NS       = 1.0;             % Interhemispheric exchange rate (years)
tau_NS_strat = 3.3;             % Interhemispheric exchange rate in the stratosphere (Fabian et al., 1968)
% Are we using the other sinks?
if ~use_other_sinks
    k_ch4_other  = 0;
    k_co_other   = 0;
    k_c2h6_other = 0;
end

%%% Convert so everything has units of days (working with julian days)
tau_NS       = tau_NS * YrToDay;
tau_NS_strat = tau_NS_strat * YrToDay;
RxNconv      = RxNconv * DaysToS;
% Which methyl chloroform reaction rate are we using?
global k_mcf_flag k_co_flag
k_mcf = k_mcf_A;
if k_mcf_flag
    k_mcf = k_mcf_B;
end
k_co = k_co_A;
if k_co_flag
    k_co = k_co_B;
end

%%% Convert some losses to consistent units 
k_ch4_other  = k_ch4_other  / YrToDay;
k_c2h6_other = k_c2h6_other / YrToDay;
k_c2h6_strat = k_c2h6_strat / YrToDay;
k_co_other   = k_co_other   / YrToDay;
k_co_strat   = k_co_strat   / YrToDay;
k_oh_strat   = k_oh_strat   / YrToDay;

%%% Guess for the initial conditions for the box model
% Troposphere
nh_12ch4 = 1575;                          % ppb
sh_12ch4 = 1520;                          % ppb
nh_13ch4 = nh_12ch4 * (1 - 47.6/1000);    % ppb
sh_13ch4 = sh_12ch4 * (1 - 47.4/1000);    % ppb
nh_mcf   = 78;                            % ppt
sh_mcf   = 68;                            % ppt
nh_n2o   = 302.0;                         % ppb
sh_n2o   = 301.0;                         % ppb
nh_c2h6  = 1200;                          % ppt
sh_c2h6  = 200;                           % ppt
nh_oh    = (gmOH/n_air)*1d9;              % ppb
sh_oh    = (gmOH/n_air)*1d9;              % ppb
nh_co    = 95;                            % ppb
sh_co    = 67;                            % ppb
% Stratosphere
nh_12ch4_S = 1570;                        % ppb
sh_12ch4_S = 1515;                        % ppb
nh_13ch4_S = nh_12ch4 * (1 - 47.5/1000);  % ppb
sh_13ch4_S = sh_12ch4 * (1 - 47.3/1000);  % ppb
nh_mcf_S   = 30;                          % ppt
sh_mcf_S   = 20;                          % ppt
nh_n2o_S   = 302.0;                       % ppb
sh_n2o_S   = 301.0;                       % ppb
nh_c2h6_S  = 300;                         % ppt
sh_c2h6_S  = 100;                         % ppt
nh_oh_S    = (1.0d6/n_air)*1d9;           % ppb
sh_oh_S    = (1.0d6/n_air)*1d9;           % ppb
nh_co_S    = 70;                          % ppb
sh_co_S    = 50;                          % ppb

% Assemble the ICs into a vector
IC = [  nh_12ch4,   sh_12ch4,   nh_13ch4,   sh_13ch4,   nh_mcf,   sh_mcf,   nh_n2o,   sh_n2o,   nh_c2h6,   sh_c2h6,   nh_oh,   sh_oh,   nh_co,   sh_co,...
      nh_12ch4_S, sh_12ch4_S, nh_13ch4_S, sh_13ch4_S, nh_mcf_S, sh_mcf_S, nh_n2o_S, sh_n2o_S, nh_c2h6_S, sh_c2h6_S, nh_oh_S, sh_oh_S, nh_co_S, sh_co_S];

%%% ODE45 parameters
Tspan = St;
opts  = odeset('MaxStep',YrToDay/12,...     % Make sure the max timestep is 1 month
               'NonNegative',1:length(IC)); % Ensure the result is positive

%%% Make a structure with the parameters
% Unit conversions
params.mm_ch4  = mm_ch4;
params.mm_mcf  = mm_mcf;
params.mm_n2o  = mm_n2o;
params.mm_c2h6 = mm_c2h6;
params.mm_oh   = mm_oh;
params.mm_co   = mm_co;
params.YrToDay = YrToDay;
params.n_air   = n_air;
params.gmOH    = gmOH;
% Rate constants/lifetimes
params.k_12ch4        = RxNconv * k_12ch4;   
params.k_13ch4        = RxNconv * k_13ch4;
params.k_mcf          = RxNconv * k_mcf;
params.k_c2h6         = RxNconv * k_c2h6;
params.k_n2o_nh       = 1/(t_n2oNH * YrToDay);
params.k_n2o_sh       = 1/(t_n2oSH * YrToDay);
params.k_ch4_strat_nh = 1/(t_ch4_strat_nh * YrToDay);
params.k_ch4_strat_sh = 1/(t_ch4_strat_sh * YrToDay);
params.k_mcf_strat    = 1/(t_mcf_strat * YrToDay);
params.k_ch4_other    = k_ch4_other;
params.k_c2h6_other   = k_c2h6_other;
params.k_c2h6_strat   = k_c2h6_strat;
params.k_co           = RxNconv * k_co;
params.k_co_strat     = k_co_strat;
params.k_co_other     = k_co_other;
params.k_oh_strat     = k_oh_strat;
params.RxNconv        = RxNconv;
params.DaysToS        = DaysToS;
params.tau_NS         = tau_NS;
params.tau_NS_strat   = tau_NS_strat;
% ODE parameters
params.Tspan   = Tspan;
params.IC      = IC;
params.odeOpts = opts;

end


%%% =======================================================================
%%% = END
%%% =======================================================================
