%%% script to tune kx 
%%% Newton Nguyen
%%% 17 July 2018

%%% Creating kx paramter vectors
kX_NH = [1.75, 2.2, 1.75, 2.2]; 
kX_SH = [2.5, 1.85, 2.5, 1.85];
% run the model 
for i = 1:length(kX_NH)
if i< 2
ignoreCO = true;
else
ignoreCO = false;
end

[emissions, concentrations] = test_kx(kX_NH(i), kX_SH(i), ignoreCO);
ems_struct{i} = emissions;
cons_struct{i} = concentrations;


mean([concentrations.nh_oh; concentrations.sh_oh])

% take ratio of both hemispheres
ratio = concentrations.nh_oh ./ concentrations.sh_oh
mean(ratio)


end
time = [1980:2016];
f = figure;
subplot(2,2, 1)
plot(time, ems_struct{1}{1}(:,1), 'r', time, ems_struct{2}{1}(:,1), 'b');
title('Methane Ems in NH Without CO')
xlabel('year')
ylabel('tg CH4')

subplot(2,2, 3)
plot(time, ems_struct{1}{1}(:,2), 'r', time, ems_struct{2}{1}(:,2), 'b');
title('Methane Ems in SH Without CO')
xlabel('year')
ylabel('tg CH4')

subplot(2,2, 2)
plot(time, ems_struct{3}{1}(:,1), 'r', time, ems_struct{4}{1}(:,1), 'b');
title('Methane Ems in NH Including CO')
xlabel('year')
ylabel('tg CH4')

subplot(2,2, 4)
plot(time, ems_struct{3}{1}(:,2), 'r', time, ems_struct{4}{1}(:,2), 'b');
title('Methane Ems in SH Including CO')
xlabel('year')
ylabel('tg CH4')
legend('NH/SH [OH]=1.2', 'NH/SH [OH]=0.8');
saveas(f, 'kX_tests.png', 'png')