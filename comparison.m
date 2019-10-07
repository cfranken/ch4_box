clear all 
close all

load('data/obs/StoredData/Turner_InputData_1980-2018_year-year.mat')
obs = out;
%c0 = ems_anal(:,1) + ems_anal(:,2);

%c0_oh = 1e6*(out_anal.nh_oh + out_anal.sh_oh)/2; 

c0_ch4 = (obs.nh_ch4 + obs.sh_ch4)/2;


c0_mcf = (obs.nh_mcf + obs.sh_mcf)/2;



%c0_nh_mcf_ems = ems(:,5);
%c0_sh_mcf_ems = ems(:,6);





load('data/obs/StoredData/Turner2017_InputData_1980-2016_year-year.mat')
obs = out;
%p0 = ems_anal(:,1) + ems_anal(:,4);
%p0 = [p0; nan; nan];
p0_ch4 = (obs.nh_ch4 + obs.sh_ch4)/2
p0_ch4 = [p0_ch4; nan; nan];
p0_mcf = (obs.nh_mcf + obs.sh_mcf)/2;
p0_mcf = [p0_mcf; nan; nan];
%p0_oh = (anal_soln{1}(:,7) + anal_soln{1}(:,8))/2*1e6;
%p0_oh = [p0_oh; nan; nan];
p0_nh_mcf_err = obs.nh_mcf_err;
p0_nh_mcf_err = [p0_nh_mcf_err; nan; nan];
p0_sh_mcf_err = obs.sh_mcf_err;
p0_sh_mcf_err = [p0_sh_mcf_err; nan; nan];
%p0_nh_mcf_ems = ems(:,3);
%p0_nh_mcf_ems = [p0_nh_mcf_ems; nan; nan];
%p0_sh_mcf_ems = ems(:,6);
%p0_sh_mcf_ems = [p0_sh_mcf_ems; nan; nan];


%%% AJT get the delta
d0_ch4 = c0_ch4 - p0_ch4
d0_mcf = c0_mcf - p0_mcf

time = [1980: 2018];

% plot figure
%figure(1)
%subplot(1,2,1)
%plot(time , c0, 'ro', time, p0, 'b-');
%title('Methane Emissions Comparison')
%xlabel('years')
%ylabel('Tg')
%legend('New Box Model', 'Old Box Model')

%subplot(122)
%plot(time, c0_oh, 'ro', time, p0_oh, 'b-')
%%hold on 
%%scatter(time, p0_oh)
%%hold on
%title('Global [OH]')
%xlabel('years')
%ylabel('molec / cm^3')

%saveas(figure(1), 'emissions_comparison.pdf', 'pdf')



figure(2)
title('Observations Comparison')
subplot(1,2,1)
plot(time, c0_ch4, 'ro', time, p0_ch4, 'b-')
%scatter(time , c0_ch4)
%hold on 
%scatter(time, p0_ch4)
ylabel('ppb')
xlabel('years')
legend('New Box Model', 'Old Box Model')

subplot(1,2,2)
plot(time, c0_mcf, 'ro', time, p0_mcf, 'b-')
%scatter(time, c0_mcf)
%hold on 
%scatter( time, p0_mcf)
title('Global [MCF]')
ylabel('ppt')
xlabel('years')
saveas(figure(2), 'obs_comparison.pdf', 'pdf')

figure(3)
title('Observations Comparison')
subplot(1,2,1)
plot(time, d0_ch4, 'ro')
ylabel('\Delta ppb')
xlabel('years')
legend('New Box Model', 'Old Box Model')

subplot(1,2,2)
plot(time, d0_mcf, 'ro')
title('Global [MCF]')
ylabel('\Delta ppt')
xlabel('years')
saveas(figure(3), 'obs_comparison_delta.pdf', 'pdf')
