%%% =======================================================================
%%% = getCO.m
%%% = Alex Turner
%%% = 04/12/2016
%%% =----------------------------------------------------------------------
%%% = NOTES
%%% =  ( 1): Reads in the carbon monoxide observations.
%%% =----------------------------------------------------------------------
%%% = INPUTS
%%% =  ( 1): dataDir -- Directory containing the data.
%%% =  ( 2): reread  -- Structure that says if we'll re-read the data.
%%% =----------------------------------------------------------------------
%%% = OUTPUTS
%%% =  ( 1): out -- A structure containing the observation information.
%%% =======================================================================

function [ out ] = getCO( dataDir, reread )

%%% Diagnostic
fprintf('   * CO\n');


%%% =======================================================================
%%% HAVE WE READ THIS DATA BEFORE?
%%% =======================================================================

%%% Build the filename
OutName = sprintf('%sobs/StoredData/co_%4i-%4i_%s-%s.mat',...
                  reread.dir,reread.sYear,reread.eYear,reread.tRes,reread.tAvg);

%%% Load pre-existing data file
if ~reread.flag
    % Check if a file exists
    if exist(OutName, 'file') == 2
        fprintf('   * LOADING OLD OBS STRUCTURE\n');
        load(OutName);
        return % Don't need to read the data
    end
end


%%% =======================================================================
%%% READ DATA
%%% =======================================================================

%%% Create the output structure
out = struct;
out.obs = struct;
out.tim = struct;
out.lat = struct;


%%% =======================================================================
%%% MOAA
%%% =======================================================================

%%% Append the directory onto the dataDir
dataDir = sprintf('%sobs/co/NOAA/month/',dataDir);

%%% Create the output structure
out = struct;
out.obs = struct;
out.tim = struct;
out.lat = struct;

%%% Define the site names, header lengths, and latitudes
% Filename structure
fNameS = 'co_%s_surface-flask_1_ccgg_month.txt';
% Make the site list with: "ls month/ch4_*_month.txt | cut -d'_' -f2 > site_list.csv"
fid = fopen(sprintf('%s/../site_list.csv',dataDir));
dat = textscan(fid,'%s %f','delimiter',',');
fclose(fid);
sNames = dat{1};
sLat   = dat{2};

%%% Read the data
for i = 1:length(sNames)
    % Current filename
    fName = sprintf('%s%s',dataDir,sprintf(fNameS,sNames{i}));
    % Get the number of header lines
    fid  = fopen(fName);
    nHDR = textscan(fid, '%s', 1,'delimiter','\n');
    nHDR = strsplit(char(nHDR{1}));
    nHDR = str2double(nHDR(end));
    fclose(fid);
    % Load the data
    dat  = importdata(fName,' ',nHDR);
    dat  = dat.data;
    tDat = datenum(dat(:,1),dat(:,2),ones(size(dat(:,1))));
    yDat = dat(:,3);
    % Remove NaNs
    ind  = ~isnan(yDat) & ~isnan(tDat);
    if sum(ind) > 0
        tDat = tDat(ind);
        yDat = yDat(ind);
        % Put the data in a structure
        out.obs.(sprintf('%s_INSTAAR',sNames{i})) = yDat;
        out.tim.(sprintf('%s_INSTAAR',sNames{i})) = tDat;
        out.lat.(sprintf('%s_INSTAAR',sNames{i})) = sLat(i);
    end
end


%%% =======================================================================
%%% SAVE THIS OBSERVATION FILE
%%% =======================================================================

%%% Save the structure
fprintf('   * SAVING OBS STRUCTURE\n');
if exist(OutName, 'file') == 2
    delete(OutName);
end
save(OutName,'out');

end


%%% =======================================================================
%%% =                             E N D                                   =
%%% =======================================================================