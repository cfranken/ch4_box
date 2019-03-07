load('case1.mat')

nh_ch4 = 1700;
sh_ch4 = 1700;
nh_co = 100;
sh_co = 100;
nh_oh = 6.6e5;
sh_oh = 6.6e5;
kx_global = 0.7;

day2sec = 60*60*24; % convert days to seconds
year2sec = 365*day2sec;
%n_air = params.n_air; % molec/cm^3 for dry atmosphere
n_air = 2.5e19;
conversion = day2sec * n_air/1d9; % conversion factor rom ppb/day to molec/cm^3 / s;
ppb2con = n_air / 1e9;


k_ch4 = 5e-15;
k_co = 2.0e-13;
%k_ch4 = params.k_12ch4; % ppb/day
%k_co =  params.k_co; % ppb/day
%k_ch4 = k_ch4 / conversion; % molec/cm^3 / s
%k_co = k_co / conversion; % molec/ cm^3 / s

kX_NH = 1;
kX_SH = 1;
%kX_NH = 1.885 ; % molec/cm^3/s
%kX_SH = 2.212 ; % molec/cm^3/s




%%% Read in concentrations 
nh_ch4 = ppb2con * nh_ch4; % molec/cm^3
sh_ch4 = ppb2con * sh_ch4; % molec/cm^3
global_ch4  = 0.5* (nh_ch4 + sh_ch4); % molec/ cm^3

nh_co = ppb2con * nh_co; % molec/ cm^3
sh_co = ppb2con * sh_co; % molec/ cm^3
global_co = 0.5* (nh_co + sh_co); % molec/ cm^3

%nh_oh = out.nh_oh; % molec/ cm^3
%sh_oh = out.sh_oh; % molec/ cm^3
global_oh = 0.5* (nh_oh + sh_oh); % molec/ cm^3


time_index = 1;
%nh_e_folds = zeros(size(out.nh_ch4));
%sh_e_folds = zeros(size(out.nh_ch4));
%global_e_folds = zeros(size(out.nh_ch4));

for t = 1:time_index

nh_jacobian = zeros(3,3);
sh_jacobian = zeros(3,3);
global_jacobian = zeros(3,3);


%%% Construct the jacobians according to Prather 1994

nh_jacobian(1,1) = -k_ch4 * nh_oh(t);
sh_jacobian(1,1) = -k_ch4 * sh_oh(t);
global_jacobian(1,1) = -k_ch4 * global_oh(t);

%%% the [2,1] element is 0


nh_jacobian(3,1) = -k_ch4 * nh_ch4(t);
sh_jacobian(3,1) = -k_ch4 * sh_ch4(t);
global_jacobian(3,1) = -k_ch4 * global_ch4(t);

% 2. Now for the CO equation 
nh_jacobian(1,2) = k_ch4 * nh_oh(t);
sh_jacobian(1,2) = k_ch4 * sh_oh(t);
global_jacobian(1,2) = k_ch4 * global_oh(t);

nh_jacobian(2,2) = -k_co * nh_oh(t);
sh_jacobian(2,2) = -k_co * sh_oh(t);
global_jacobian(2,2) = -k_co * global_oh(t);

nh_jacobian(3,2) = k_ch4 * nh_ch4(t) - k_co * nh_co(t);
sh_jacobian(3,2) = k_ch4 * sh_ch4(t) - k_co * sh_co(t);
global_jacobian(3,2) = k_ch4 * global_ch4(t) - k_co * global_co(t);

% 3. The Oh reaactions, d OH / dt 
nh_jacobian(1,3) = -k_ch4 * nh_oh(t);
sh_jacobian(1,3) = -k_ch4 * sh_oh(t);
global_jacobian(1,3) = -k_ch4 * global_oh(t);

nh_jacobian(2,3) = -k_co * nh_oh(t);
sh_jacobian(2,3) = -k_co * sh_oh(t);
global_jacobian(2,3) = -k_co * global_oh(t);

nh_jacobian(3,3) = -k_ch4 * nh_ch4(t) - k_co * nh_co(t) - kX_NH;
sh_jacobian(3,3) = -k_ch4 * sh_ch4(t) - k_co * sh_co(t) - kX_SH;
global_jacobian(3,3) = -k_ch4 * global_ch4(t) - k_co * global_co(t) - kx_global;

nh_jacobian = nh_jacobian';
sh_jacobian = sh_jacobian';
global_jacobian = global_jacobian';


%%% Take the eigenvalues 
D_nh = eig(nh_jacobian);
D_sh = eig(sh_jacobian);
D_global = eig(global_jacobian);



% CH4 lifetimes
result.ch4_nh_lifetime(t) = -1/D_nh(2)/year2sec;
result.ch4_sh_lifetime(t) = -1/D_sh(2)/year2sec;
result.ch4_global_lifetime(t) = -1/D_global(2)/year2sec;
result.ch4_ss(t) = -1 ./(global_jacobian(1,1)*year2sec);

% CO lifetimes
result.co_nh_lifetime(t) = -1/D_nh(3)/day2sec;
result.co_sh_lifetime(t) = -1/D_sh(3)/day2sec;
result.co_global_lifetime(t) = -1/D_global(3)/day2sec;
result.co_ss(t) = -1 ./(global_jacobian(2,2)*day2sec);

% OH lifetimes 
result.oh_nh_lifetime(t) = -1/D_nh(1);
result.oh_sh_lifetime(t) = -1/D_sh(1);
result.oh_global_lifetime(t) = -1/D_global(1);
result.oh_ss(t) = -1 ./global_jacobian(3,3);


end
