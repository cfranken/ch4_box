ignoreCO= true;
kX_NH = 1.75;
kX_SH = 2.5;
[emissions, concentrations] = test_kx(kX_NH, kX_SH,ignoreCO);



mean([concentrations.nh_oh; concentrations.sh_oh])

% take ratio of both hemispheres
ratio = concentrations.nh_oh ./ concentrations.sh_oh;

mean(ratio)