function obs_struct = useMopit(St, obs_struct, tRes)

disp('*** Replacing surface data with MOPIT obs***')

mopit_start = [2000, 1, 1];
St_years = datevec(St); St_years = St_years(:,1);
ind_start = find( St_years==mopit_start); % the starting index for MOPIT data 


sYear = 2000; % beginning of MOPIT record
eYear = min( datenum(2017, 1, 1), St(end));
eYear = datevec(eYear); eYear = eYear(1);

 % end of MOPiT record
tRes = 'month';
St_mopit    = getTime(sYear,eYear,tRes); % Time vector

tRes = 'year';
St_blockOutput    = getTime(sYear,eYear,tRes); % Time vector

mopit = xlsread('mopit_co.xlsx');

nh_co = mopit(:,4);
nh_co_err = mopit(:,3);
fDays = 365.25;

sh_co = mopit(:,7);
sh_co_err = mopit(:,6);

% Take annual averages of the data 
nh_co = annualAvg(nh_co);
sh_co = annualAvg(sh_co);
nh_co_err = annualAvg(nh_co_err);
sh_co_err = annualAvg(sh_co_err);




% replace the data
% Also rescale with surface data
obs_struct.nh_co(ind_start : ind_start + length(nh_co_err)-1) = 0.9377*nh_co;
obs_struct.sh_co(ind_start:ind_start+length(sh_co_err)-1) = 0.8574*sh_co;
obs_struct.nh_co_err(ind_start : ind_start + length(nh_co_err)-1) = 0.1224*nh_co_err;
obs_struct.sh_co_err(ind_start:ind_start+length(sh_co_err)-1) = 0.1469*sh_co_err;



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

