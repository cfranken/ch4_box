%%% =======================================================================
%%% = getCOems.m
%%% = Alex Turner
%%% = 05/02/2017
%%% =----------------------------------------------------------------------
%%% = NOTES
%%% =  ( 1): Create the CO emissions and puts them onto our temporal 
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

function [ out ] = getOHems( St, tRes, dataDir )

%%% Diagnostic
fprintf('   * OH\n');

%%% Which emissions do we want to use?
const_ems  = true;
lin_change = false;
if const_ems
    % Total OH source (constant)
    tot_oh  = 4500; % Tg/yr
    frac_nh = 0.50; % Fraction of anthropogenic emissions in the NH
    % OH Emissions
    ems_nh = zeros(size(St)) + tot_oh * frac_nh;
    ems_sh = zeros(size(St)) + tot_oh * (1 - frac_nh);
end
if lin_change
    frac_nh = 0.50; % Fraction of emissions in the NH
    date_S  = datenum(1980,1,1);
    date_E  = datenum(2010,1,1);
    ems_S   = 9000;    % Tg/yr in 1980
    ems_E   = 9500;    % Tg/yr in 2010
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
