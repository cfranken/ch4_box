function obs_struct = useMopit(St, obs_struct, tRes)
%%% Takes surface obs and replaces with MOPIT obs
%%% Rescales surface data, prior to MOPIT, according to MOPIT scaling
%%% inputs:
%%% 1. St: Our time vector
%%% 2. obs_struct: The concentration obs structure from makeObs.m
%%% 3. tRes: time resolution. Either 'year' or 'month'
%%% output:
%%%obs_struct: the processed observation struct with replaced MOPIT data and rescaled surface data 
%%% Newton Nguyen 4 July 2019


disp('*** Replacing surface data with MOPIT obs***')

mopit_start = [2000, 1, 1];
St_years = datevec(St); St_years = St_years(:,1);
ind_start = find( St_years==mopit_start); % the starting index for MOPIT data 


sYear = 2000; % beginning of MOPIT record
eYear = min( datenum(2017, 1, 1), St(end));
eYear = datevec(eYear); eYear = eYear(1);

% Read the MOPIT datafile
mopit = xlsread('mopit_co.xlsx');

% assign data to names 
nh_co = mopit(:,4);
nh_co_err = mopit(:,3);


sh_co = mopit(:,7);
sh_co_err = mopit(:,6);

% Take annual averages of the data 
nh_co = annualAvg(nh_co);
sh_co = annualAvg(sh_co);
nh_co_err = annualAvg(nh_co_err);
sh_co_err = annualAvg(sh_co_err);

% compute scaling factors by only using timeseries during MOPIT era
nh_co_scale = mean(nh_co ./ obs_struct.nh_co(ind_start : ind_start + length(nh_co_err)-1))
sh_co_scale = mean(sh_co ./ obs_struct.sh_co(ind_start:ind_start+length(sh_co_err)-1))
nh_co_err_scale = mean(nh_co_err ./ obs_struct.nh_co_err(ind_start : ind_start + length(nh_co_err)-1))
sh_co_err_scale = mean(sh_co_err ./ obs_struct.sh_co_err(ind_start:ind_start+length(sh_co_err)-1))

% rescale the pre-MOPIT surface obs with MOPIT obs
obs_struct.nh_co(1:ind_start - 1) = nh_co_scale * obs_struct.nh_co(1:ind_start - 1);
obs_struct.sh_co(1 : ind_start-1) = sh_co_scale * obs_struct.sh_co(1 : ind_start-1);
obs_struct.nh_co_err(1 : ind_start - 1) = nh_co_err_scale * obs_struct.nh_co_err(1 : ind_start - 1);
obs_struct.sh_co_err(1 : ind_start - 1) = sh_co_err_scale * obs_struct.sh_co_err(1 : ind_start - 1);

% Replace Surface record corresponding to MOPIT times with MOPIT obs
obs_struct.nh_co(ind_start : ind_start + length(nh_co_err)-1) = nh_co;
obs_struct.sh_co(ind_start:ind_start+length(sh_co_err)-1) = sh_co;
obs_struct.nh_co_err(ind_start : ind_start + length(nh_co_err)-1) = nh_co_err;
obs_struct.sh_co_err(ind_start:ind_start+length(sh_co_err)-1) = sh_co_err;

%%%% End of function

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

