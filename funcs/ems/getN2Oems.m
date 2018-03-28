%%% =======================================================================
%%% = getN2Oems.m
%%% = Alex Turner
%%% = 05/02/2017
%%% =----------------------------------------------------------------------
%%% = NOTES
%%% =  ( 1): Create the N2O emissions and puts them onto our temporal 
%%% =        grid.
%%% =----------------------------------------------------------------------
%%% = INPUTS
%%% =  ( 1): St      -- Our time vector.
%%% =  ( 2): tRes    -- String containing the temporal resolution.
%%% =  ( 3): dataDir -- Directory containing the data.
%%% =----------------------------------------------------------------------
%%% = OUTPUTS
%%% =  ( 1): out -- Structure containing the methane emissions.
%%% =======================================================================

function [ out ] = getN2Oems( St, tRes, dataDir )

%%% Diagnostic
fprintf('   * N2O\n');

%%% Which emissions do we want to use?
edgar        = false;
const_ems    = false;
lin_increase = true;
if edgar
    % Read in EDGAR
    fname       = sprintf('%s/ems/n2o/%s',dataDir,'EDGARv42FT2012_N2O.xls');
    [dat, ~, ~] = xlsread(fname,'N2O_timeseries');
    tDat        = dat(8,:);     % Year
    tDat        = datenum(tDat,ones(size(tDat)),ones(size(tDat)));
    yDat        = dat(238,:);   % Global emissions (Gg/yr)
    % Interpolate it to my temporal grid
    if strcmp(tRes,'year') || strcmp(tRes,'YEAR') || strcmp(tRes,'yearly')
        fDays        = 365; % Number of days in the block average
        [tDat, yDat] = BlockAverage(tDat,yDat,ones(size(tDat)),fDays);
    end
    ind  = ~isnan(tDat) & ~isnan(yDat);
    oDat = interp1(tDat(ind),yDat(ind),St,'spline');
    oDat(isnan(oDat)) = nanmax(oDat);
    % Convert it to the units we need and add a natural source
    date_S       = datenum(1980,1,1);
    date_E       = datenum(2010,1,1);
    nat_source_S = 10000;      % Gg/yr in 1980
    nat_source_E = 14000;      % Gg/yr in 2010
    nat_source_T = (nat_source_E - nat_source_S)/(date_E - date_S);
    nat_source_O = nat_source_S - (nat_source_T*date_S);
    nat_source   = nat_source_T*St + nat_source_O; % Tg/yr
    frac_nh_anth = 0.9;     % Fraction of anthropogenic emissions in the NH
    frac_nh_nat  = 0.4;     % Fraction of natural emissions in the NH
    ems_nh       = oDat*frac_nh_anth     + nat_source*frac_nh_nat;
    ems_sh       = oDat*(1-frac_nh_anth) + nat_source*(1-frac_nh_nat);
    ems_nh       = ems_nh * 1d-3;   % Tg/yr
    ems_sh       = ems_sh * 1d-3;   % Tg/yr
end
if const_ems
    % Total N2O source (constant)
    tot_n2o = 10;   % Tg/yr
    frac_nh = 0.75; % Fraction of emissions in the NH
    % N2O Emissions
    ems_nh = zeros(size(St)) + tot_n2o * frac_nh;
    ems_sh = zeros(size(St)) + tot_n2o * (1 - frac_nh);
end
if lin_increase
    frac_nh = 0.75; % Fraction of emissions in the NH
    date_S  = datenum(1980,1,1);
    date_E  = datenum(2010,1,1);
    ems_S   = 10;      % Tg/yr in 1980
    ems_E   = 14;      % Tg/yr in 2010
    ems_T   = (ems_E - ems_S)/(date_E - date_S);
    ems_O   = ems_S - (ems_T*date_S);
    ems     = ems_T*St + ems_O; % Tg/yr
    ems_nh  = ems * frac_nh;
    ems_sh  = ems - ems_nh;
end
% Make the structure
out.nh = ems_nh;
out.sh = ems_sh;

end


%%% =======================================================================
%%% =                             E N D                                   =
%%% =======================================================================