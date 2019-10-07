if use_MOPIT_CO
sYear = 2000; % beginning of MOPIT record
eYear = min( datenum(2017, 1, 1), St(1));
 % end of MOPiT record
tRes = 'month';
St_mopit    = getTime(sYear,eYear,tRes); % Time vector

tRes = 'year';
St_blockOutput    = getTime(sYear,eYear,tRes); % Time vector

mopit = xlsread('mopit_co.xlsx')

nh_co = mopit(:,4);
nh_co_err = mopit(:,3);

sh_co = mopit(:,7);
sh_co_err = mopit(:,6);
    [tDat, yDat] = BlockAverage_CO(St_blockOutput,St_mopit,nh_co,ones(size(St_mopit)),fDays);
    [tDat, nh_coerr] = BlockAverage_CO(St_blockOutput,St_mopit,nh_co_err,ones(size(St_mopit)),fDays);
    [tDat, sh_co] = BlockAverage_CO(St_blockOutput,St_mopit,sh_co,ones(size(St_mopit)),fDays);
    [tDat, sh_coerr] = BlockAverage_CO(St_blockOutput,St_mopit,sh_co_err,ones(size(St_mopit)),fDays);


end
fDays=365.25;
    [tDat, yDat] = BlockAverage(St,nh_co,ones(size(St)),fDays);


% 
% let's separate the MOPIT co with the rest of the CO data. The makeobs will read all the field names, so let's not make the confusing 