%%% a script to plot impacts of fixed OH concentration v. source in our test

% load the data files 
load('case1_test.mat')
load('case2_test.mat')
time = [1980:2016];
case1.ch4_ems = case1.ems{1}(:,1) + case1.ems{1}(:,2);
case2.ch4_ems = case2.ems{1}(:,1) + case2.ems{1}(:,2);

%%% OH concentrations 

case1.nh_oh_con = case1.concentrations.nh_oh;
case2.nh_oh_con = case2.concentrations.nh_oh;


%%% OH source
case1.oh_ems= case1.ems{1}(:,11) + case1.ems{1}(:,12);

case2.oh_ems= case2.ems{1}(:,11) + case1.ems{1}(:,12);



figure
subplot(2,2,1)
plot(time, case1.ch4_ems, time, case2.ch4_ems);
xlabel('year');
ylabel('tg CH4');
title('Methane Emissions');

subplot(2,2,2);
plot(time, case1.nh_oh_con', time, case2.nh_oh_con);
xlabel('year')
ylabel('molecules/cm^3')
title('NH OH concentration')

subplot(2,2,3);
plot(time, case1.oh_ems, time, case2.oh_ems)
xlabel('year')
ylabel('tg OH')
title('OH source')

legend('case 1', 'Case 2')
