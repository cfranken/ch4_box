%%% =======================================================================
%%% = readData.m
%%% = Alex Turner
%%% = 06/06/2016
%%% =----------------------------------------------------------------------
%%% = NOTES
%%% =  ( 1): Reads a csv file that we made.
%%% =----------------------------------------------------------------------
%%% = INPUTS
%%% =  ( 1): baseName -- Prefix for the plots.
%%% =----------------------------------------------------------------------
%%% = OUTPUTS
%%% =   N/A
%%% =  ( 1): St       -- Our time vector.
%%% =  ( 2): obs      -- Structure containing the observations.
%%% =  ( 3): ems      -- Emission sources (and OH) for the box model.
%%% =  ( 4): mod      -- Structure containing the model results.
%%% =  ( 5): IC       -- Initial conditions.
%%% =======================================================================

function [ St, obs, mod, ems, IC ] = readData( baseName )

%%% Read the data
dat = importdata(sprintf(baseName,'Data'),',',2);
dat = dat.data;
% Fill the years
St = dat(:,1);
% Fill obs
obs.nh_ch4    = dat(:,2);
obs.sh_ch4    = dat(:,3);
obs.nh_ch4c13 = dat(:,4);
obs.sh_ch4c13 = dat(:,5);
obs.nh_mcf    = dat(:,6);
obs.sh_mcf    = dat(:,7);
obs.nh_n2o    = dat(:,8);
obs.sh_n2o    = dat(:,9);
obs.nh_c2h6   = dat(:,10);
obs.sh_c2h6   = dat(:,11);
obs.nh_co     = dat(:,12);
obs.sh_co     = dat(:,13);
% Fill obs error
obs.nh_ch4_err    = dat(:,14);
obs.sh_ch4_err    = dat(:,15);
obs.nh_ch4c13_err = dat(:,16);
obs.sh_ch4c13_err = dat(:,17);
obs.nh_mcf_err    = dat(:,18);
obs.sh_mcf_err    = dat(:,19);
obs.nh_n2o_err    = dat(:,20);
obs.sh_n2o_err    = dat(:,21);
obs.nh_c2h6_err   = dat(:,22);
obs.sh_c2h6_err   = dat(:,23);
obs.nh_co_err     = dat(:,24);
obs.sh_co_err     = dat(:,25);
% Fill model output
mod.nh_ch4    = dat(:,26);
mod.sh_ch4    = dat(:,27);
mod.nh_ch4c13 = dat(:,28);
mod.sh_ch4c13 = dat(:,29);
mod.nh_mcf    = dat(:,30);
mod.sh_mcf    = dat(:,31);
mod.nh_n2o    = dat(:,32);
mod.sh_n2o    = dat(:,33);
mod.nh_c2h6   = dat(:,34);
mod.sh_c2h6   = dat(:,35);
mod.nh_oh     = dat(:,36);
mod.sh_oh     = dat(:,37);
mod.nh_co     = dat(:,38);
mod.sh_co     = dat(:,39);
% Fill emissions
ems.nh_ch4    = dat(:,40);
ems.sh_ch4    = dat(:,41);
ems.nh_ch4c13 = dat(:,42);
ems.sh_ch4c13 = dat(:,43);
ems.nh_mcf    = dat(:,44);
ems.sh_mcf    = dat(:,45);
ems.nh_n2o    = dat(:,46);
ems.sh_n2o    = dat(:,47);
ems.nh_c2h6   = dat(:,48);
ems.sh_c2h6   = dat(:,49);
ems.nh_oh     = dat(:,50);
ems.sh_oh     = dat(:,51);
ems.nh_co     = dat(:,52);
ems.sh_co     = dat(:,53);
ems.tau_TS    = dat(:,54);  % strat-trop exchange
ems.kX_NH     = dat(:,55);  % OH loss rate
ems.kX_SH     = dat(:,56);  % OH loss rate

%%% Read the ICs
dat = importdata( sprintf(baseName,'ICs'),',',2);
dat = dat.data;
IC  = dat;

end


%%% =======================================================================
%%% = END
%%% =======================================================================
