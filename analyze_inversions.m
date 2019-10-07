clf 
close all
clear all

% Newton Nguyen 
% 05/14/2018
% A script to analyze the data from inversions % Exporting results of test to test_results.mat for python plotting because pythong plots are being used in the paper 

% Reference Case from AJT 2017 paper
load('case0_test_tau0')

case0.prior = ems;
case0.concentrations = out_anal;
case0.ems = anal_soln{1};
case0.nh_ch4_ems = anal_soln{1}(:,1);
case0.sh_ch4_ems = anal_soln{1}(:,2);
case0.nh_co_ems = anal_soln{1}(:,13);
case0.sh_co_ems = anal_soln{1}(:,14);
case0.nh_oh = anal_soln{1}(:,11);
case0.sh_oh = anal_soln{1}(:,12);
ems_struct = makeEmsStruct(ems_anal);
prior = makeEmsStruct(ems);
save('case0_test_tau0')

%save('case0_test.mat')

%%%1. Run the basecase
%%% Alex's PNAS Paper 
% ignoreCO = true
% interactiveOH = false
% fixedOH = true 
load('case1_test_tau0')


% Give variables names 
case1.concentrations = out_anal;
case1.ems = anal_soln{1};
case1.nh_ch4_ems = anal_soln{1}(:,1);
case1.sh_ch4_ems = anal_soln{1}(:,2);
case1.nh_co_ems = anal_soln{1}(:,13);
case1.sh_co_ems = anal_soln{1}(:,14);
case1.nh_oh_ems = anal_soln{1}(:,11);
case1.sh_oh_ems = anal_soln{1}(:,12);
ems_struct = makeEmsStruct(ems_anal);
prior = makeEmsStruct(ems);

delta_nh_ch4 = out_anal.nh_ch4 - obs.nh_ch4;
delta_sh_ch4 = out_anal.sh_ch4 - obs.sh_ch4;
time = [1980 : 1980 + length(St)-1];

%figure(1)
%subplot(2,2,1)
%plot(time , obs.nh_ch4, 'go', time, out_anal.nh_ch4, 'k-')
%title('NH CH4 concentrations')
%
%subplot(2,2,3)
%plot(time , obs.sh_ch4, 'go', time, out_anal.sh_ch4, 'k-')
%title('SH CH4 concentrations')
%legend('Observations', 'Model Output')

%subplot(2,2,2)
%plot(time, delta_nh_ch4, 'ro')

%title('Model - Obs in NH')

%subplot(2,2,4)
%plot(time, delta_nh_ch4, 'bo')
%title('Model - Obs in SH')
%saveas(figure(1), 'ch4_fit.pdf', 'pdf')




save('case1_test_tau0')


load('case2_test_tau0')


% Give variables names 
case2.concentrations = out_anal;
case2.ems = anal_soln;
case2.nh_ch4_ems = anal_soln{1}(:,1);
case2.sh_ch4_ems = anal_soln{1}(:,2);
case2.nh_co_ems = anal_soln{1}(:,13);
case2.sh_co_ems = anal_soln{1}(:,14);
case2.nh_oh_ems = anal_soln{1}(:,11);
case2.sh_oh_ems = anal_soln{1}(:,12);
ems_struct = makeEmsStruct(ems_anal);
prior = makeEmsStruct(ems);
save('case2_test_tau0')


%%% Case 3:
% ignoreCO = false
% interactiveOH = true
% fixedOH = true


load('case3_test_tau0')


% Give variables names 
case3.concentrations = out_anal;
case3.ems = anal_soln;
case3.nh_ch4_ems = anal_soln{1}(:,1);
case3.sh_ch4_ems = anal_soln{1}(:,2);
case3.nh_co_ems = anal_soln{1}(:,13);
case3.sh_co_ems = anal_soln{1}(:,14);
case3.nh_oh_ems = anal_soln{1}(:,11);
case3.sh_oh_ems = anal_soln{1}(:,12);
ems_struct = makeEmsStruct(ems_anal);
prior = makeEmsStruct(ems);

save('case3_test_tau0')


load('case4_test_tau0')

case4.concentrations = out_anal;
case4.ems = anal_soln;
case4.nh_ch4_ems = anal_soln{1}(:,1);
case4.sh_ch4_ems = anal_soln{1}(:,2);
case4.nh_co_ems = anal_soln{1}(:,13);
case4.sh_co_ems = anal_soln{1}(:,14);
case4.nh_oh_ems = anal_soln{1}(:,11);
case4.sh_oh_ems = anal_soln{1}(:,12);
ems_struct = makeEmsStruct(ems_anal);
prior = makeEmsStruct(ems);

save('case4_test_tau0')

load('case5_test_tau0')

case5.concentrations = out_anal;
case5.ems = anal_soln;
case5.nh_ch4_ems = anal_soln{1}(:,1);
case5.sh_ch4_ems = anal_soln{1}(:,2);
case5.nh_co_ems = anal_soln{1}(:,13);
case5.sh_co_ems = anal_soln{1}(:,14);
case5.nh_oh_ems = anal_soln{1}(:,11);
case5.sh_oh_ems = anal_soln{1}(:,12);
ems_struct = makeEmsStruct(ems_anal);
prior = makeEmsStruct(ems);
save('case5_test_tau0')

case0.ch4 = case0.nh_ch4_ems + case0.sh_ch4_ems;
case1.ch4 = case1.nh_ch4_ems + case1.sh_ch4_ems;
case2.ch4 = case2.nh_ch4_ems + case2.sh_ch4_ems;
case3.ch4 = case3.nh_ch4_ems + case3.sh_ch4_ems;
case4.ch4 = case4.nh_ch4_ems + case4.sh_ch4_ems;
case5.ch4 = case5.nh_ch4_ems + case5.sh_ch4_ems;

fprintf('case 0 CH4 ems at %4.0f in year 2000 and %4.0f in year 2017 \n \n', case0.ch4(21), case0.ch4(end))
fprintf('case 1 CH4 ems at %4.0f in year 2000 and %4.0f in year 2017 \n', case1.ch4(21), case1.ch4(end))
fprintf('case 2 CH4 ems at %4.0f in year 2000 and %4.0f in year 2017 \n', case2.ch4(21), case2.ch4(end))
fprintf('case 3 CH4 ems at %4.0f in year 2000 and %4.0f in year 2017 \n', case3.ch4(21), case3.ch4(end))
fprintf('case 4 CH4 ems at %4.0f in year 2000 and %4.0f in year 2017 \n', case4.ch4(21), case4.ch4(end))
fprintf('case 5 CH4 ems at %4.0f in year 2000 and %4.0f in year 2017 \n', case5.ch4(21), case5.ch4(end))

figure(2)
plot(time, case0.ems(:,5), 'ro', time, case1.ems(:,5), 'go', time, ems(:,5), 'ko')
legend