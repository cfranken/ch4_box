%%% =======================================================================
%%% = getC2H6ems.m
%%% = Alex Turner
%%% = 05/02/2017
%%% =----------------------------------------------------------------------
%%% = NOTES
%%% =  ( 1): Create the methane emissions and puts them onto our temporal 
%%% =        grid.  Can either use constant emissions, constant with a 30
%%% =        Tg/yr jump in 2007, or EDGAR v4.2FT2010.
%%% =----------------------------------------------------------------------
%%% = INPUTS
%%% =  ( 1): St      -- Our time vector.
%%% =  ( 2): tRes    -- String containing the temporal resolution.
%%% =  ( 3): dataDir -- Directory containing the data.
%%% =----------------------------------------------------------------------
%%% = OUTPUTS
%%% =  ( 1): out -- Structure containing the methane emissions.
%%% =======================================================================

function [ out ] = getC2H6ems( St, tRes, dataDir )

%%% Diagnostic
fprintf('   * C2H6\n');

%%% Which emissions do we want to use?
const_ems = true;
if const_ems
    % Total C2H6 source (constant)
    ems_globe = 13000; % Gg/yr (Xiao et al)
    ems_ratio = 0.8; % 80% of emissions are in the NH
    ems_nh    = ems_globe * ems_ratio; % Gg/yr
    ems_sh    = ems_globe - ems_nh;
    % C2H6 Emissions
    ems_nh = zeros(size(St)) + ems_nh;
    ems_sh = zeros(size(St)) + ems_sh;
end
% Make the structure
out.nh = ems_nh;
out.sh = ems_sh;

end


%%% =======================================================================
%%% =                             E N D                                   =
%%% =======================================================================