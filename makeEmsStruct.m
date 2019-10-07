function ems_struct = makeEmsStruct(ems_array)

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

ems_struct = struct;

ems_struct.nh_ch4 = ems_array(:,1);
ems_struct.sh_ch4 = ems_array(:,2);
ems_struct.nh_ch4c13 = ems_array(:,3);
ems_struct.sh_ch4c13 = ems_array(:,4);
ems_struct.nh_mcf = ems_array(:,5);
ems_struct.sh_mcf = ems_array(:,6);
ems_struct.nh_n2o = ems_array(:,7);
ems_struct.sh_n2o = ems_array(:,8);
ems_struct.nh_c2h6 = ems_array(:,9);
ems_struct.sh_c2h6 = ems_array(:,10);
ems_struct.nh_oh = ems_array(:,11);
ems_struct.sh_oh = ems_array(:,12);
ems_struct.nh_co = ems_array(:,13);
ems_struct.sh_co = ems_array(:,14);
ems_struct.tau = ems_array(:,15);
ems_struct.kx_nh = ems_array(:,16);
ems_struct.kx_sh = ems_array(:,17);



