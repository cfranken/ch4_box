% Newton Nguyen 
% 05/14/2018
% A script to run tests of the inversions located in DriverScript
% Exporting results of test to test_results.mat for python plotting because pythong plots % are prettier


%%%1. Run the basecase
%%% Alex's PNAS Paper 
% ignoreCO = true
% interactiveOH = false
% fixedOH = true 
DriverScript_case1


% Give variables names 
case1.concentrations = out_anal;
case1.ems = anal_soln;
case1.nh_ch4_ems = anal_soln{1}(:,1);
case1.sh_ch4_ems = anal_soln{1}(:,2);
case1.nh_co_ems = anal_soln{1}(:,13);
case1.sh_co_ems = anal_soln{1}(:,14);
case1.nh_oh_ems = anal_soln{1}(:,11);
case1.sh_oh_ems = anal_soln{1}(:,12);

save('case1_test.mat')
DriverScript_case2


% Give variables names 
case2.concentrations = out_anal;
case2.ems = anal_soln;
case2.nh_ch4_ems = anal_soln{1}(:,1);
case2.sh_ch4_ems = anal_soln{1}(:,2);
case2.nh_co_ems = anal_soln{1}(:,13);
case2.sh_co_ems = anal_soln{1}(:,14);
case2.nh_oh_ems = anal_soln{1}(:,11);
case2.sh_oh_ems = anal_soln{1}(:,12);
save('case2_test.mat')


%%% Case 3:
% ignoreCO = false
% interactiveOH = true
% fixedOH = true


DriverScript_case3


% Give variables names 
case3.concentrations = out_anal;
case3.ems = anal_soln;
case3.nh_ch4_ems = anal_soln{1}(:,1);
case3.sh_ch4_ems = anal_soln{1}(:,2);
case3.nh_co_ems = anal_soln{1}(:,13);
case3.sh_co_ems = anal_soln{1}(:,14);
case3.nh_oh_ems = anal_soln{1}(:,11);
case3.sh_oh_ems = anal_soln{1}(:,12);

save('case3_test.mat')

DriverScript_case4;

case4.concentrations = out_anal;
case4.ems = anal_soln;
case4.nh_ch4_ems = anal_soln{1}(:,1);
case4.sh_ch4_ems = anal_soln{1}(:,2);
case4.nh_co_ems = anal_soln{1}(:,13);
case4.sh_co_ems = anal_soln{1}(:,14);
case4.nh_oh_ems = anal_soln{1}(:,11);
case4.sh_oh_ems = anal_soln{1}(:,12);

save('case4_test.mat')

DriverScript_case5

case5.concentrations = out_anal;
case5.ems = anal_soln;
case5.nh_ch4_ems = anal_soln{1}(:,1);
case5.sh_ch4_ems = anal_soln{1}(:,2);
case5.nh_co_ems = anal_soln{1}(:,13);
case5.sh_co_ems = anal_soln{1}(:,14);
case5.nh_oh_ems = anal_soln{1}(:,11);
case5.sh_oh_ems = anal_soln{1}(:,12);

save('case5_test.mat')


