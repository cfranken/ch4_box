%%% =======================================================================
%%% = writeData.m
%%% = Alex Turner
%%% = 06/06/2016
%%% =----------------------------------------------------------------------
%%% = NOTES
%%% =  ( 1): Saves the output as a csv file.
%%% =----------------------------------------------------------------------
%%% = INPUTS
%%% =  ( 1): St       -- Our time vector.
%%% =  ( 2): obs      -- Structure containing the observations.
%%% =  ( 3): ems      -- Emission sources (and OH) for the box model.
%%% =  ( 4): mod      -- Structure containing the model results.
%%% =  ( 5): IC       -- Initial conditions.
%%% =  ( 6): baseName -- Prefix for the plots.
%%% =----------------------------------------------------------------------
%%% = OUTPUTS
%%% =   N/A
%%% =======================================================================

function [ ] = writeData( St, obs, mod, ems, IC, baseName )

%%% Convert the matrix to a structure
ems = disassembleEms(ems);

%%% Assemble all the data into a single matrix
% Initialized
nObs = 10;          % NH/SH CH4, NH/SH d13C, NH/SH MCF, NH/SH N2O, NH/SH C2H6
nEms = 12;          % NH/SH CH4, NH/SH d13C, NH/SH MCF, NH/SH N2O, NH/SH C2H6, NH/SH OH
nSt  = length(St);
dat  = nan(nSt,nObs+nObs+nObs+nEms+1);
% Fill the years
dat(:,1)  = St;
% Fill obs
dat(:,2)  = obs.nh_ch4;
dat(:,3)  = obs.sh_ch4;
dat(:,4)  = obs.nh_ch4c13;
dat(:,5)  = obs.sh_ch4c13;
dat(:,6)  = obs.nh_mcf;
dat(:,7)  = obs.sh_mcf;
dat(:,8)  = obs.nh_n2o;
dat(:,9)  = obs.sh_n2o;
dat(:,10) = obs.nh_c2h6;
dat(:,11) = obs.sh_c2h6;
dat(:,12) = obs.nh_co;
dat(:,13) = obs.sh_co;
% Fill obs error
dat(:,14) = obs.nh_ch4_err;
dat(:,15) = obs.sh_ch4_err;
dat(:,16) = obs.nh_ch4c13_err;
dat(:,17) = obs.sh_ch4c13_err;
dat(:,18) = obs.nh_mcf_err;
dat(:,19) = obs.sh_mcf_err;
dat(:,20) = obs.nh_n2o_err;
dat(:,21) = obs.sh_n2o_err;
dat(:,22) = obs.nh_c2h6_err;
dat(:,23) = obs.sh_c2h6_err;
dat(:,24) = obs.nh_co_err;
dat(:,25) = obs.sh_co_err;
% Fill model output
dat(:,26) = mod.nh_ch4;
dat(:,27) = mod.sh_ch4;
dat(:,28) = mod.nh_ch4c13;
dat(:,29) = mod.sh_ch4c13;
dat(:,30) = mod.nh_mcf;
dat(:,31) = mod.sh_mcf;
dat(:,32) = mod.nh_n2o;
dat(:,33) = mod.sh_n2o;
dat(:,34) = mod.nh_c2h6;
dat(:,35) = mod.sh_c2h6;
dat(:,36) = mod.nh_oh;
dat(:,37) = mod.sh_oh;
dat(:,38) = mod.nh_co;
dat(:,39) = mod.sh_co;
% Fill emissions
dat(:,40) = ems.nh_ch4;           % NH CH4
dat(:,41) = ems.sh_ch4;           % SH CH4
dat(:,42) = ems.nh_ch4c13;        % NH CH4C13
dat(:,43) = ems.sh_ch4c13;        % SH CH4C13
dat(:,44) = ems.nh_mcf;           % NH MCF
dat(:,45) = ems.sh_mcf;           % SH MCF
dat(:,46) = ems.nh_n2o;           % NH N2O
dat(:,47) = ems.sh_n2o;           % SH N2O
dat(:,48) = ems.nh_c2h6;          % NH C2H6
dat(:,49) = ems.sh_c2h6;          % SH C2H6
dat(:,50) = ems.nh_oh;            % NH OH
dat(:,51) = ems.sh_oh;            % SH OH
dat(:,52) = ems.nh_co;            % NH CO
dat(:,53) = ems.sh_co;            % SH CO
dat(:,54) = ems.tau_TS;           % strat-trop exchange
dat(:,55) = ems.kX_NH;            % NH OH arbitrary loss rate
dat(:,56) = ems.kX_SH;            % SH OH arbitrary loss rate
% Make the strings
dat_Head = {'JULIAN DATE',...
            'NH CH4 (obs)','SH CH4 (obs)','NH d13C (obs)','SH d13C (obs)','NH MCF (obs)','SH MCF (obs)','NH N2O (obs)','SH N2O (obs)','NH C2H6 (obs)','SH C2H6 (obs)','NH CO (obs)','SH CO (obs)',...
            'NH CH4 (err)','SH CH4 (err)','NH d13C (err)','SH d13C (err)','NH MCF (err)','SH MCF (err)','NH N2O (err)','SH N2O (err)','NH C2H6 (err)','SH C2H6 (err)','NH CO (err)','SH CO (err)',...
            'NH CH4 (mod)','SH CH4 (mod)','NH d13C (mod)','SH d13C (mod)','NH MCF (mod)','SH MCF (mod)','NH N2O (mod)','SH N2O (mod)','NH C2H6 (mod)','SH C2H6 (mod)','NH OH (mod)','SH OH (mod)','NH CO (mod)','SH CO (mod)',...
            'NH CH4 (ems)','SH CH4 (ems)','NH d13C (ems)','NH d13C (ems)','NH MCF (ems)','SH MCF (ems)','NH N2O (ems)','SH N2O (ems)','NH C2H6 (ems)','SH C2H6 (ems)','NH OH (ems)','SH OH (ems)','NH CO (ems)','SH CO (ems)','tau_TS (ems)','NH kX (ems)','SH kX (ems)'};
dat_Unit = {'days since Jan-1-0000',...
                     'ppb',         'ppb',       'permil',       'permil',         'ppt',         'ppt',         'ppb',         'ppb',          'ppt',          'ppt',        'ppb',        'ppb',...
                     'ppb',         'ppb',       'permil',       'permil',         'ppt',         'ppt',         'ppb',         'ppb',          'ppt',          'ppt',        'ppb',        'ppb',...
                     'ppb',         'ppb',       'permil',       'permil',         'ppt',         'ppt',         'ppb',         'ppb',          'ppt',          'ppt',  'molec/cm3',  'molec/cm3',        'ppb',        'ppb',...
                   'Tg/yr',       'Tg/yr',       'permil',       'permil',       'Gg/yr',       'Gg/yr',       'Tg/yr',       'Tg/yr',        'Gg/yr',        'Gg/yr',      'Tg/yr',      'Tg/yr',      'Tg/yr',      'Tg/yr',          'yr',       's^-1',       's^-1'};
fspecA = '%s,';
fspecB = '%9i,';
for i = 2:size(dat,2)
    fspecA = sprintf('%s%%s,',fspecA);
    fspecB = sprintf('%s%%9.3f,',fspecB);
end
fspecA = sprintf('%s\\n',fspecA);
fspecB = sprintf('%s\\n',fspecB);

%%% Write the data
fid = fopen(sprintf(baseName,'Data'),'w');
fprintf(fid,fspecA,dat_Head{:}); % Header
fprintf(fid,fspecA,dat_Unit{:}); % Units
% Data
for i = 1:nSt
    fprintf(fid,fspecB,dat(i,:));
end
fclose(fid);

%%% Write the ICs
fid = fopen(sprintf(baseName,'ICs'),'w');
fprintf(fid,'NH 12CH4,SH 12CH4,NH 13CH4,SH 13CH4,NH MCF,SH MCF,NH N2O,SH N2O,NH C2H6,SH C2H6,NH OH,SH OH,NH CO,SH CO,NH-S 12CH4,SH-S 12CH4,NH-S 13CH4,SH-S 13CH4,NH-S MCF,SH-S MCF,NH-S N2O,SH-S N2O,NH-S C2H6,SH-S C2H6,NH-S OH,SH-S OH,NH-S CO,SH-S CO,\n');  % Header
fprintf(fid,'ppb,ppb,permil,permil,ppt,ppt,ppb,ppb,ppt,ppt,ppb,ppb,ppb,ppb,ppb,ppb,permil,permil,ppt,ppt,ppb,ppb,ppt,ppt,ppb,ppb,ppb,ppb,\n');                                            % Units
fprintf(fid,'%9.3f,%9.3f,%9.3f,%9.3f,%9.3f,%9.3f,%9.3f,%9.3f,%9.3f,%9.3f,%9.3f,%9.3f,%9.3f,%9.3f,%9.3f,%9.3f,%9.3f,%9.3f,%9.3f,%9.3f,%9.3f,%9.3f,%9.3f,%9.3f,%9.3f,%9.3f,%9.3f,%9.3f,\n',IC);                   % Data
fclose(fid);

end


%%% =======================================================================
%%% = END
%%% =======================================================================
