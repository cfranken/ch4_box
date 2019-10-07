%%% =======================================================================
%%% = DriverScript.m
%%% = Alex Turner
%%% = Originally created on 04/12/2016
%%% =----------------------------------------------------------------------
%%% = NOTES:
%%% =
%%% = This is the driver script for the 2-box model methane inversion.
%%% = There are currently two different inversions implemented: (1) a
%%% = linear or non-linear deterministic inversion following Rodgers (2000)
%%% = and (2) an inversion using the non-linear Covariance Matrix Adaptation
%%% = Evolution Strategy (CMA-ES).  Case (1) requires us to compute gradients
%%% = and only allows for Gaussian errors.  Case (2) is  a stochastic method
%%% = that automatically tunes the proposal distribution to improve the sampling,
%%% = however it does not provide error statistics that are consistent with
%%% = the distributions.  Case (2) also allows one to specify non-analytic
%%% = distributions (e.g., bounded Gaussians or uniform distributions).
%%% =======================================================================


%%
%%% =======================================================================
%%% 1. Initialize
%%% =======================================================================
profile off
%%% Clear the MatLab space
clf
clear all
close all
clc

%%% Header
fprintf('\n ***********************************\n')
fprintf(' *** STARTING GLOBAL 2-BOX MODEL ***\n')
fprintf(' ***********************************\n')

%%% Define the directories
baseDir = pwd;
utilDir = sprintf('%s/funcs/', baseDir);
dataDir = sprintf('%s/data/',  baseDir);
outDir  = sprintf('%s/output/',baseDir);

%%% Add the utility functions
addpath(utilDir);
addpath(sprintf('%s/obs',               utilDir));
addpath(sprintf('%s/ems',               utilDir));
addpath(sprintf('%s/model',             utilDir));
addpath(sprintf('%s/util',              utilDir));
addpath(sprintf('%s/plot',              utilDir));
addpath(sprintf('%s/inv',               utilDir));
addpath(sprintf('%s/inv/deterministic', utilDir));
addpath(sprintf('%s/inv/stochastic',    utilDir));

%%% Define the time period
sYear = 1980;
eYear = 2017;
%eYear = 2100;
tRes  = 'year';     % Can be 'year' or 'month' (year preferred)
tAvg  = 'year';     % Smooth the observations
St    = getTime(sYear,eYear,tRes); % Time vector
nT    = length(St);

%%% Export variables to mat file
export_data = true; % do we want to export data to data_filename.mat?
data_filename  = 'MOPIT_test';

%%% Describing experiment to be exported to .mat file 
disp('*** Commencing the CO test with MOPIT versus Surface CO ***')

%%% Execute in parallel?
run_parallel = false;
if run_parallel
    nWorkers     = 4;
    setupParallel(run_parallel,nWorkers);
end


%%% What kind of inversions do we want to do?
do_deterministic = true;    % Rodgers (2000)
do_cmaes         = false;    % Covariance Matrix Adaptation Evolution Strategy

%%% For reading the observations
% Do we want to reread the raw data?
reread.flag  = true;
% Other flags for re-reading
reread.sYear = sYear;
reread.eYear = eYear;
reread.tRes  = tRes;
reread.tAvg  = tAvg;
reread.dir   = dataDir;

%%% Other options and flags
% Use globals for some flags
global fixedCH4 fixedOH onlyCH4 onlyMCF schaefer          % Linear inversion
global k_mcf_flag smooth_MCF set_MCF_EMS MCF_EMS_val      % Methyl Chloroform
global k_co_flag use_strat interactive_OH use_other_sinks ignoreCO % Other
global use_MOPIT_CO
global no_temporal_correlation large_prior % inversion tests on prior constraints
% Plotting flags
ftype           = 'pdf';    % Type of plots to make? (eps, pdf, tif, or png)
plot_prior      = false;     % Plot the prior?
plot_raw        = false;    % Plot the raw observations?
plot_old_cmaes  = false;    % Plot an old CMA-ES solution (false means run a new one)
% General flags
use_strat       = false;     % Use a stratosphere?
interactive_OH  = true;     % Allow OH feedbacks?
use_other_sinks = false;     % Use non-OH sinks?
% Linear inversion flags
use_other_sinks = false;     % Use non-OH sinks?
% Linear inversion flags
det_linear      = false;     % Use a linear deterministic inversion?
fixedCH4        = false;    % Use fixed methane emissions
fixedOH         = true;    % Use fixed OH anomalies
onlyCH4         = false;    % Only invert for methane emissions
ignoreCO = false; % keep CO emissions fixed
onlyMCF         = false;    % Only invert for MCF emissions
schaefer        = false;    % Case that is most similar to Schaefer et al.
% Flags for priors in inversions
no_temporal_correlation = false; % Run with no temporal correlation? Should be run with large_prior
large_prior = false; % Run with large prior in emissions? 

% MCF sensitivity test flags
k_co_flag       = true;     % Use k_CO that AJT derived
k_mcf_flag      = true;     % Use k_MCF that AJT derived
smooth_MCF      = false;    % Smooth the MCF emissions with a 5-year filter?
set_MCF_EMS     = false;    % Set post-2000 emissions to a fixed value?
MCF_EMS_val     = 0.0;      % Fixed post-2000 MCF emissions value (Gg/yr)
reduce_MCFerr   = false;    % Reduce the errors in MCF observations?
MCF_ERR_val     = 2.0;      % Error in MCF observations (ppt)
% Flags for other tests to run
use_OH_stratMLO = false;    % Use the OH derived from MLO strat ozone?
use_Ed          = false;    % Use Ed Dlugokencky's hemispheric averages?
use_Turner_Bootstrap = false; % use data from Turner et al, 2017?
use_MOPIT_CO = true; % Use MOPIT data?

%%% Set the seed for repeatability
rng('default');


%%
%%% =======================================================================
%%% 2. Load the obs
%%% =======================================================================

%%% Diagnostic
fprintf('\n *** LOADING THE OBSERVATIONS *** \n');

%%% Load the observations
% Structures with with three fields:
% - "obs":  Observations from each NOAA site (ppb)
% - "tim":  Julian date for the observation
% - "lat":  Latitude of the NOAA site
try % Add a try-catch statement in case the user hasn't downloaded the data
    ch4_obs     = getCH4(dataDir,reread);      % CH4 observations (ppb)
    ch4c13_obs  = getCH4C13(dataDir,reread);   % delta13C observations (permil)
    ch4h2_obs   = getCH4H2(dataDir,reread);    % deltaD observations (permil)
    mcf_obs     = getMCF(dataDir,reread);      % Methylchloroform observations (ppt)
    n2o_obs     = getN2O(dataDir,reread);      % N2O observations (ppb)
    c2h6_obs    = getC2H6(dataDir,reread);     % Ethane observations (ppt)
    co_obs      = getCO(dataDir,reread);       % carbon monoxide observations (ppb)
    o3strat_obs = getO3strat(dataDir,reread);  % Stratospheric ozone observations (DU)
catch % Some data is missing
    try % See if ethane is the only problem
        fprintf(' * SOME DATA IS MISSING\n');
        ch4_obs     = getCH4(dataDir,reread);
        ch4c13_obs  = getCH4C13(dataDir,reread);
        ch4h2_obs   = getCH4H2(dataDir,reread);
        mcf_obs     = getMCF(dataDir,reread);
        n2o_obs     = getN2O(dataDir,reread);
        co_obs      = getCO(dataDir,reread);
        o3strat_obs = getO3strat(dataDir,reread);
        c2h6_obs    = NaN;
    catch % Otherwise, set the observation structures to NaN
        fprintf(' * UNABLE TO READ OBSERVATIONS!\n');
        ch4_obs     = NaN;
        ch4c13_obs  = NaN;
        ch4h2_obs   = NaN;
        mcf_obs     = NaN;
        n2o_obs     = NaN;
        c2h6_obs    = NaN;
        co_obs      = NaN;
        o3strat_obs = NaN;
    end
end

%%% Make the observation structure
% Structure with 12 fields:
% - NH/SH CH4    obs & err (ppb)
% - NH/SH CH4C13 obs & err (permil)
% - NH/SH MCF    obs & err (ppt)
% - NH/SH N2O    obs & err (ppb)
% - NH/SH C2H6   obs & err (ppt)
% - NH/SH CO     obs & err (ppb)
obs_mopit = makeObs(St,tAvg,ch4_obs,ch4c13_obs,mcf_obs,n2o_obs,c2h6_obs,co_obs,dataDir,reread);

if use_MOPIT_CO
obs_mopit = useMopit(St, obs_mopit, tRes);
end

%

if use_Turner_Bootstrap
    turnerFname = sprintf('%sobs/StoredData/Turner_InputData_%4i-%4i_%s-%s.mat',...
                  dataDir,reread.sYear,reread.eYear,reread.tRes,reread.tAvg);
    ajt_obs = load(turnerFname);
    obs     = ajt_obs.out;

% Get rid of CO data before 1990
coYear = datenum(1991, 1, 1);
ind = find(St<coYear);
obs.nh_co(ind(1) : ind(end)) = nan;
obs.sh_co(ind(1) : ind(end)) = nan;

end


% blow up CO error:
%obs.nh_co_err(:)=500;
%obs.sh_co_err(:)=500;
%%% Use Ed Dlugokencky's obs? (sensitivity test)
if use_Ed
    ajt_obs = obs;
    ed_obs  = getEdObs(dataDir,ajt_obs,St,tAvg,reread);
    obs     = ed_obs;
    plotEdObs(St,ajt_obs,ed_obs,sprintf('%s/%s/raw_EdObs.%s',outDir,tRes,ftype))
end

%%% Reduce MCF errors?  (sensitivity test)
if reduce_MCFerr
    obs.nh_mcf_err = min([obs.nh_mcf_err,MCF_ERR_val*ones(size(obs.nh_mcf_err))],[],2);
    obs.sh_mcf_err = min([obs.sh_mcf_err,MCF_ERR_val*ones(size(obs.sh_mcf_err))],[],2);
end

%%% Diagnostics (check the raw data)
if plot_raw
    deseasonalize  = true;
    plot_all_sites = false;
    plotAllObs(St,obs,ch4_obs,   tAvg, 'ch4' ,sprintf('%s/%s/raw_%%s_%%s.%s',outDir,tRes,ftype),deseasonalize,plot_all_sites);
    plotAllObs(St,obs,ch4c13_obs,tAvg, 'd13C',sprintf('%s/%s/raw_%%s_%%s.%s',outDir,tRes,ftype),deseasonalize,plot_all_sites);
    %plotAllObs(St,obs,ch4h2_obs, tAvg, 'dD',  sprintf('%s/%s/raw_%%s_%%s.%s',outDir,tRes,ftype),deseasonalize,plot_all_sites);
    plotAllObs(St,obs,mcf_obs,   tAvg, 'mcf' ,sprintf('%s/%s/raw_%%s_%%s.%s',outDir,tRes,ftype),deseasonalize,plot_all_sites);
    plotAllObs(St,obs,n2o_obs,   tAvg, 'n2o' ,sprintf('%s/%s/raw_%%s_%%s.%s',outDir,tRes,ftype),deseasonalize,plot_all_sites);
    plotAllObs(St,obs,c2h6_obs,  tAvg, 'c2h6',sprintf('%s/%s/raw_%%s_%%s.%s',outDir,tRes,ftype),deseasonalize,plot_all_sites);
    plotAllObs(St,obs,co_obs,    tAvg, 'co'  ,sprintf('%s/%s/raw_%%s_%%s.%s',outDir,tRes,ftype),deseasonalize,plot_all_sites);
end


%%
%%% =======================================================================
%%% 3. Load the emissions (all will be arrays with a length matching "St")
%%% =======================================================================

%%% Diagnostic
fprintf('\n *** LOADING THE EMISSIONS *** \n');

%%% Get the CH4 emissions
% Stucture with two fields
% - "nh": CH4 emissions from the Northern hemisphere (Tg/yr)
% - "sh": CH4 emissions from the Southern hemisphere (Tg/yr)
ch4_ems = getCH4ems(St,tRes,dataDir);

%%% Get the delta13C composition for NH/SH CH4 emissions
% Stucture with two fields
% - "nh": delta13C composition from the Northern hemisphere (permil)
% - "sh": delta13C composition from the Southern hemisphere (permil)
ch4c13_ems = getCH4C13ems(St,tRes,dataDir);

%%% Get the MCF emissions (assumed to be in NH only)
% Stucture with two fields
% - "prinn":     MCF emissions from Prinn (Gg/yr)
% - "mcculloch": MCF emissions from McCulloch (Gg/yr)
mcf_ems = getMCFems(St,tRes,dataDir);

%%% Get the N2O emissions
% Stucture with two fields
% - "nh": N2O emissions from the Northern hemisphere (Tg/yr)
% - "sh": N2O emissions from the Southern hemisphere (Tg/yr)
n2o_ems = getN2Oems(St,tRes,dataDir);

%%% Get the C2H6 emissions
% Stucture with two fields
% - "nh": C2H6 emissions from the Northern hemisphere (Tg/yr)
% - "sh": C2H6 emissions from the Southern hemisphere (Tg/yr)
c2h6_ems = getC2H6ems(St,tRes,dataDir);

%%% Get the OH emissions
% Stucture with two fields
% - "nh": OH emissions from the Northern hemisphere (Tg/yr)
% - "sh": OH emissions from the Southern hemisphere (Tg/yr)
oh_ems = getOHems(St,tRes,dataDir);

%%% Get the CO emissions
% Stucture with two fields
% - "nh": CO emissions from the Northern hemisphere (Tg/yr)
% - "sh": CO emissions from the Southern hemisphere (Tg/yr)
co_ems = getCOems(St,tRes,dataDir);
% somehow, the NH emissions are way too low:
%co_ems.nh(:) = 1400;

%%
%%% =======================================================================
%%% 4. Initialize the 2-box model
%%% =======================================================================

%%% Diagnostic
fprintf('\n *** RUN THE 2-BOX MODEL WITH PRIOR FLUXES *** \n');

%%% OH scaling factor
oh_scale.nh = ones(nT,1);
oh_scale.sh = ones(nT,1);
% Derive OH from the stratospheric ozone?
if use_OH_stratMLO
    OH_sensitivity  = 4.2;          % a 1% increase in strat O3 leads to a 4.2% decrease in OH (Murray et al., 2013)
    O3_site         = 'mlo_NOAA';   % which site to use?
    fDays           = 365.25*2;     % How long of a smoothing?
    tO3             = o3strat_obs.tim.(O3_site);
    yO3             = o3strat_obs.obs.(O3_site);
    yO3             = DeseasonalizeData(tO3,yO3,fDays);
    [tO3, yO3, ~]   = BlockAverage_AltError(tO3,yO3,ones(size(tO3)),365.25);
    oh_change       = yO3 / nanmean(yO3); % Convert strat O3 to OH change
    oh_change       = 1 ./ ((oh_change - 1) * OH_sensitivity + 1);
    yOH             = interp1(tO3,oh_change,St);
    yOH(isnan(yOH)) = 1;
    % Store this OH
    oh_scale.nh = yOH;
    oh_scale.sh = yOH;
end

%%% Strat-trop exchange
tau_TS = 9.0 * ones(nT,1); % years
if ~use_strat
    % Set this to something high, Inf results in trouble:
    %tau_TS(:) = Inf; % No exchange with stratosphere
    tau_TS(:) = 1e4;
end

%%% Arbitrary reactions with OH
% CF Needed to adapt NH as there would otherwise be a rather large IH
% difference in OH
kX_NH = 0.99*ones(nT,1); % s^-1 for 6600 tg/yr OH source
kX_SH = 1.23*ones(nT,1); % s^-1

%%% Structure of sources with 17 fields:
% - NH CH4 emissions
% - SH CH4 emissions
% - NH CH4C13 composition
% - SH CH4C13 composition
% - NH MCF emissions
% - SH MCF emissions
% - NH N2O emissions
% - SH N2O emissions
% - NH C2H6 emissions
% - SH C2H6 emissions
% - NH OH emissions
% - SH OH emissions
% - NH CO emissions
% - SH CO emissions
% - Strat-trop exchange
% - NH arbitrary OH reaction rate
% - SH arbitrary OH reaction rate
ems.nh_ch4    = ch4_ems.nh;
ems.sh_ch4    = ch4_ems.sh;
ems.nh_ch4c13 = ch4c13_ems.nh;
ems.sh_ch4c13 = ch4c13_ems.sh;
ems.nh_mcf    = mcf_ems.nh;
ems.sh_mcf    = mcf_ems.sh;
ems.nh_n2o    = n2o_ems.nh;
ems.sh_n2o    = n2o_ems.sh;
ems.nh_c2h6   = c2h6_ems.nh;
ems.sh_c2h6   = c2h6_ems.sh;
ems.nh_oh     = oh_ems.nh;
ems.sh_oh     = oh_ems.sh;
ems.nh_co     = co_ems.nh;
ems.sh_co     = co_ems.sh;
ems.tau_TS    = tau_TS;
ems.kX_NH     = kX_NH;
ems.kX_SH     = kX_SH;
% Convert the structure to a matrix
ems = assembleEms(ems);

%%% Run the box model
params = getParameters(St); % Only need to do this once
IC     = params.IC;         % Guess for the inital conditions
out    = boxModel_wrapper(St,ems,IC,params);
if plot_prior
    plotNewObs(St,out,obs,sprintf('%s/%s/prior_%%s.%s',outDir,tRes,ftype));
    %writeData(St,obs,out,ems,IC,sprintf('%s/%s/prior_%%s.csv',outDir,tRes));
    %plotObs(St,out,obs,sprintf('%s/%s/prior_%%s.%s',outDir,tRes,ftype));
    %plotDrivers(St,ems,NaN*ems,sprintf('%s/%s/prior_%%s.%s',outDir,tRes,ftype),dataDir);
end


%%
%%% =======================================================================
%%% 5. Deterministic inversion (Rodgers, 2000)
%%% =======================================================================

if do_deterministic
    
    %%% Diagnostic
    fprintf('\n *** DETERMINISTIC INVERSION *** \n');
    
    %%% Invert
    [anal_soln,jacobian_ems,jacobian_IC,reltol,abstol, mati] = invert_methane(St,obs_mopit,ems,IC,params,det_linear,run_parallel);


% amkes obs struct with surface observations 
obs_surface = makeObs(St,tAvg,ch4_obs,ch4c13_obs,mcf_obs,n2o_obs,c2h6_obs,co_obs,dataDir,reread);

    [anal_soln2,jacobian_ems2,jacobian_IC2,reltol2,abstol2, mati2] = invert_methane(St,obs_surface,ems,IC,params,det_linear,run_parallel);


    
    
    %%% Plot the Jacobians
    %[jacobian_ems,jacobian_IC] = define_Jacobian( St, ems, IC, params, run_parallel );


%    plotJacobian(St,jacobian_ems,tRes,sprintf('%s/%s/jacobian_%%s.%s',outDir,tRes,ftype));
    
    %%% Try plotting the solution
    ems_anal = anal_soln{1};
    IC_anal  = anal_soln{2};
    % Comment this out for now:
    out_anal = boxModel_wrapper(St,ems_anal,IC_anal,params);
%    plotNewObs(St,out_anal,obs,sprintf('%s/%s/anal_%%s.%s',outDir,tRes,ftype));
    %writeData(St,obs,out_anal,ems_anal,IC_anal,sprintf('%s/%s/anal_%%s.csv',outDir,tRes));
    %plotObs(St,out_anal,obs,sprintf('%s/%s/anal_%%s.%s',outDir,tRes,ftype));
    %plotDrivers(St,ems_anal,ems,sprintf('%s/%s/anal_%%s.%s',outDir,tRes,ftype),dataDir);
    
end

% plot the original mopit obs 

mopit = xlsread('mopit_co.xlsx');

nh_co = mopit(:,4);
nh_co_err = mopit(:,3);


sh_co = mopit(:,7);
sh_co_err = mopit(:,6);

% Take annual averages of the data 
nh_co = annualAvg(nh_co);
sh_co = annualAvg(sh_co);
nh_co_err = annualAvg(nh_co_err);
sh_co_err = annualAvg(sh_co_err);

if export_data
    fprintf('Exporting all variables in this run to %s \n', data_filename)
    save(data_filename);
end

%%% Plot the results 
clf
close all


figure(8)
time = [1: length(St)];

subplot(211)
plot(time, obs_surface.nh_co, 'ro', time, obs_mopit.nh_co, 'bo')
title('NH CO obs rescaled')
xlabel('years')
ylabel('ppb')
legend('Surface Obs', 'MOPIT Obs');

subplot(212)
plot(time, obs_surface.sh_co, 'ro', time, obs_mopit.sh_co, 'bo')
title('SH CO obs rescaled')

saveas(figure(8), 'mopit_co_scaled.png', 'png')
figure(7)
subplot(211)
plot(time, anal_soln{1}(:,13), 'ro', time, anal_soln2{1}(:,13), 'bo')
title('NH CH4 Ems')
ylabel('Tg')
xlabel('years')

subplot(212)
plot(time, anal_soln{1}(:,14), 'ro', time, anal_soln2{1}(:,14), 'bo')
title('SH CH4 Ems')
ylabel('Tg')
xlabel('years')
legend('MOPIT CO', 'Surface CO')
saveas(figure(7), 'mopit_co_ems', 'png')


% plot the original mopit obs 

mopit = xlsread('mopit_co.xlsx');

nh_co = mopit(:,4);
nh_co_err = mopit(:,3);


sh_co = mopit(:,7);
sh_co_err = mopit(:,6);

% Take annual averages of the data 
nh_co = annualAvg(nh_co);
sh_co = annualAvg(sh_co);
nh_co_err = annualAvg(nh_co_err);
sh_co_err = annualAvg(sh_co_err);
time = [2000:2017];


figure(9)
subplot (211)
plot(time, obs_surface.nh_co(21:end), 'ro', time , nh_co, 'go')
xlabel('years')
ylabel('ppb')
title('Unscaled NH CO obs')

subplot(212)
plot(time, obs_surface.sh_co(21:end), 'ro', time , sh_co, 'go');
xlabel('years')
ylabel('ppb')
title('Unscaled SH CO obs')
saveas(figure(9), 'mopit_co_unscaled', 'png')

function annual_avg = annualAvg(data)
%%% Takes the annual average of a timeseries 
%%% inputs:
%%% 1. data: The (nX1) array to take an average over
%%% outputs:
%%% 1. mean value (double)

n_blocks = floor(length(data)/12);
annual_avg = nan(n_blocks,1);
for i = 1:n_blocks
ind1 = 12*(i-1) + 1;
ind2 = ind1 + 11;
annual_avg(i) = nanmean(data(ind1:ind2));
end
end
