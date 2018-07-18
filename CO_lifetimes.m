%%% Plotting CO lifetimes 
%%% Newton Nguyen
%%% July 7, 2018

case5_out = calculate_CO_lifetime(case5, 'case5_test.mat');



function out = calculate_CO_lifetime(case_name, datafile)
% case_name = name of case 
% datafile: name of .mat file e.g. filename.mat (str)
k_co = 2e-13;
sec2day = 24*3600; % convert seconds to days

load(datafile);
nh_data = case_name.concentrations.nh_oh;
sh_data = case_name.concentrations.sh_oh;
;
out.nh_lifetime = (1/k_co * sec2day).*(1./nh_data);
out.sh_lifetime = (1/k_co*sec2day) .* (1./sh_data);
end

