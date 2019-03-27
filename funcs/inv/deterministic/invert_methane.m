%%% = Alex Turner
%%% = 04/12/2016
%%% =----------------------------------------------------------------------
%%% = NOTES
%%% =  ( 1): Performs a deterministic inversion.
%%% =  ( 2): Depending on the flags, it can do a linear inversion or a
%%% =        non-linear inversion.  The non-linear inversion can use either
%%% =        Gauss-Newton or Levenberg-Marquardt.
%%% =  ( 3): Prior covariance matrix is defined in the subfunction called:
%%% =        "update_solution".cd ,,.
%%% =----------------------------------------------------------------------
%%% = INPUTS
%%% =  ( 1): St           -- Our time vector.
%%% =  ( 2): obs          -- Structure with the observations.
%%% =  ( 3): ems_p        -- Matrix with the prior emissions (and OH).
%%% =  ( 4): IC_p         -- Vector with the prior initial conditions.
%%% =  ( 5): params       -- Structure with parameters for the box model.
%%% =  ( 6): linear       -- Are we doing a linear inversion?  (True/False)
%%% =  ( 7): run_parallel -- Are we running in parallel?  (True/False).
%%% =----------------------------------------------------------------------
%%% = OUTPUTS
%%% =  ( 1): soln    -- Cell array with the emissions and ICs.
%%% =  ( 2): K_ems   -- Jacobian for the emission sources (and OH).
%%% =  ( 3): K_IC    -- Jacobian for the ICs.
%%% =  ( 4): reldiff -- Vector of relative differences while iterating.
%%% =  ( 5): absdiff -- Vector of absolute differences while iterating.
%%% =======================================================================

function [ soln, K_ems, K_IC, reldiff, absdiff, matr] = invert_methane( St, obs, ems_p, IC_p, params, linear, run_parallel )

%%% Define tolerances
reltol  = 1e-6;
abstol  = 1e-6;
kmax    = 5;
k       = 1;
iter    = true;
reldiff = nan(kmax,1);
absdiff = nan(kmax,1);

%%% Initialize the Levenberg-Marquardt parameters (gamma = zero reverts to IMAP)
LM_param.chi2  = NaN;
LM_param.gamma = 10;
LM_param.gamma = 1;

%%% Are we doing a linear or nonlinear inversion?
if linear
    LM_param.gamma = 0;
    iter           = false;
end

%%% Diagnostics
if linear
    fprintf('   * STARTING LINEAR INVERSION\n');
else
    if LM_param.gamma == 0
        fprintf('   * STARTING NONLINEAR IMAP INVERSION\n');
    else
        fprintf('   * STARTING NONLINEAR LEVENBERG-MARQUARDT INVERSION\n');
    end
end

%%% Get the solution for the first step
[K_ems,K_IC]    = define_Jacobian(St,ems_p,IC_p,params,run_parallel);
[soln,LM_param, matr] = update_solution(St,ems_p,IC_p,ems_p,IC_p,LM_param,K_ems,K_IC,obs,params);
disp('so far so good')
while iter
    % Update solution
    ems_i = soln{1};
    IC_i  = soln{2};
%    ems_i(24,:)
%    ems_p(24,:)
%    plot(IC_i-IC_p)

    chi2o = LM_param.chi2;
    % Get new solution
    [K_ems,K_IC]    = define_Jacobian(St,ems_i,IC_i,params,run_parallel);
    
    [soln,LM_param, matr] = update_solution(St,ems_i,IC_i,ems_p,IC_p,LM_param,K_ems,K_IC,obs,params);
    % save iteration steps
    FF(:,k) = matr.F;
    % Get new relative difference
    absdiff_ems = abs(soln{1}(:) - ems_i(:));
    absdiff_IC  = abs(soln{2}(:) -  IC_i(:));
    reldiff_ems = absdiff_ems ./ min(abs(soln{1}(:)),abs(ems_i(:)));
    reldiff_IC  = absdiff_IC  ./ min(abs(soln{2}(:)),abs( IC_i(:)));
    absdiff(k)  = sum([absdiff_ems;absdiff_IC]);
    reldiff(k)  = rms([reldiff_ems;reldiff_IC]);
    % Check convergence
    if (kmax <= k)
        iter = false;
    %elseif (10*LM_param.chi2 < (numel(ems_p)+numel(IC_p)))
    %    iter = false;
    elseif LM_param.chi2 ~= chi2o % Did we accept the previous solution?
        if (reldiff(k) < reltol) || (absdiff(k) < abstol)
            iter = false;
        end
    end
    fprintf('ITERATION %3i/%3i: reldiff = %3.2e & absdiff = %3.2e\n',k,kmax,reldiff(k),absdiff(k))
    k = k + 1;
end
%matr.FF = FF;
end

function [ out, LM_param, matr ] = update_solution( St, ems_i, IC_i, ems_p, IC_p, LM_param, jacobian_ems, jacobian_IC, obs, params )

%%% Alternate cases to run
global fixedCH4 fixedOH onlyCH4 onlyMCF schaefer use_strat ignoreMCF ignoreCO
global fitKX
if onlyCH4
    obs.nh_ch4c13(:) = NaN;
    obs.sh_ch4c13(:) = NaN;
    obs.nh_mcf(:)    = NaN;
    obs.sh_mcf(:)    = NaN;
    obs.nh_co(:)     = NaN;
    obs.sh_co(:)     = NaN;
end
if onlyMCF
    obs.nh_ch4c13(:) = NaN;
    obs.sh_ch4c13(:) = NaN;
    obs.nh_ch4(:)    = NaN;
    obs.sh_ch4(:)    = NaN;
    obs.nh_co(:)     = NaN;
    obs.sh_co(:)     = NaN;
end
if schaefer
    obs.nh_mcf(:) = NaN;
    obs.sh_mcf(:) = NaN;
    obs.nh_co(:)  = NaN;
    obs.sh_co(:)  = NaN;
    fixedOH       = true;
end
% Don't use ethane observations for inversion
obs.nh_c2h6(:) = NaN;
obs.sh_c2h6(:) = NaN;

%%% Get dimensions
nY = size(jacobian_ems,1);
nT = length(St);
nE = size(ems_i,2);
nI = size(jacobian_IC,2);

%%% Assemble the state and observation vectors
F  = assembleObs(boxModel_wrapper(St,ems_i,IC_i,params));
y  = assembleObs(obs);
xa = assembleStateVector(ems_p,IC_p);
xi = assembleStateVector(ems_i,IC_i);

%%% Determine which rows to throw out (no observations!)
indGood = ~isnan(y);
y       = y(indGood);
F       = F(indGood);

%%% Assemble the jacobians
% Sources
K_ems = nan(nY,nT*nE);
%if ignoreCO
%    jacobian_ems(:,:,13:14)=0;
%end


ii    = 1;
for i = 1:nE
    for j = 1:nT
        K_ems(:,ii) = jacobian_ems(:,j,i);
        ii          = ii + 1;
    end
end
size(K_ems)
% ICs
K_IC = jacobian_IC; % This one is already assembled
% Full jacobian
K = [K_ems,K_IC];
% Remove NaNs
K = K(indGood,:);

%%% Create the prior error covariance matrix
% Components
Sa_ch4     =    20^2*ones(nT,1);
Sa_ch4c13  =    10^2*ones(nT,1);
%Sa_mcf_nh  = max([.2*ems_p(:,5),1.5*ones(nT,1)],[],2).^2;
Sa_mcf_nh  = max([.02*ems_p(:,5),0.15*ones(nT,1)],[],2).^2;
Sa_mcf_sh  =   0.15^2*ones(nT,1);
Sa_n2o     =   2.0^2*ones(nT,1);
Sa_c2h6    =  5000^2*ones(nT,1);
Sa_oh      =  250^2*ones(nT,1);
Sa_co      =   300^2*ones(nT,1);
Sa_tau     =   3.0^2*ones(nT,1); 
Sa_IC      =    [30,30,10,10,15,15,5,5,100,100,...
                 30,30,10,10,15,15,5,5,100,100].^2;
% CF: Let's just go lazy here and use 5% of the IC as prior uncertainty for
% now:
Sa_IC = (eps*params.IC).^2;
tau_ch4    = 5; % yr
tau_ch4c13 = 5; % yr
tau_mcf    = 3; % yr
tau_n2o    = 5; % yr
tau_c2h6   = 1; % yr
tau_oh     = 3; % yr
tau_co     = 3; % yr
tau_tau    = 0; % yr
% Alternate cases
if fixedCH4
Sa_ch4 = eps^2*ones(nT,1); % Fixed CH4
end
if fixedOH
    Sa_oh = eps^2*ones(nT,1); % Fixed OH
end
if onlyCH4
    Sa_oh     = eps^2*ones(nT,1); % Fixed OH
    Sa_ch4c13 = eps^2*ones(nT,1); % Fixed d13C
    Sa_mcf_nh = eps^2*ones(nT,1); Sa_mcf_sh = Sa_mcf_nh; % Fixed MCF
end
if onlyMCF
    Sa_ch4    = eps^2*ones(nT,1); % Fixed CH4
    Sa_ch4c13 = eps^2*ones(nT,1); % Fixed d13C
end
if schaefer
    Sa_oh     = eps^2*ones(nT,1); % Fixed OH
    Sa_mcf_nh = eps^2*ones(nT,1); Sa_mcf_sh = Sa_mcf_nh; % Fixed MCF
end
% CF: Just create some arbitrary prior uncertainty for the arbitrary
% reaction rates here as well, as those were missing:
Sa_kx_NH = eps^2*ones(nT,1); % Fixed K
Sa_kx_SH = eps^2*ones(nT,1); % Fixed K 
% We can add flags later if we want to fit kx instead of OH sources!
if fitKX
    disp('Fitting KX')
    Sa_kx_NH = 0.92^2*ones(nT,1); 
    Sa_kx_SH = 0.92^2*ones(nT,1); 
end

% Construct matrix
if ignoreCO
    obs.nh_co_err = obs.nh_co_err*1e20;
    obs.sh_co_err = obs.sh_co_err*1e20;
    Sa_co = eps^2*ones(nT,1);
end

Sa_ems = [Sa_ch4,Sa_ch4,Sa_ch4c13,Sa_ch4c13,Sa_mcf_nh,Sa_mcf_sh,...
          Sa_n2o,Sa_n2o,Sa_c2h6,  Sa_c2h6,  Sa_oh,    Sa_oh,...
	  Sa_co, 0.07*Sa_co, Sa_tau, Sa_kx_NH, Sa_kx_SH];
Sa     = diag(assembleStateVector(Sa_ems,Sa_IC));


tau    = [tau_ch4,tau_ch4,tau_ch4c13,tau_ch4c13,tau_mcf,tau_mcf,...
          tau_n2o,tau_n2o,tau_c2h6,  tau_c2h6,  tau_oh, tau_oh, tau_co, tau_co, tau_tau]*365.25;
Sa     = fillDiagonalsAnal(Sa,tau,St);


%%% Construct the observational error covariance matrix
So = [obs.nh_ch4_err;obs.sh_ch4_err;obs.nh_ch4c13_err;obs.sh_ch4c13_err;obs.nh_mcf_err;obs.sh_mcf_err;...
      obs.nh_n2o_err;obs.sh_n2o_err;obs.nh_c2h6_err;  obs.sh_c2h6_err;obs.nh_co_err;obs.sh_co_err].^2;
So = diag(So(indGood));

%%% Get inverse covariance matrices
if any(tau > 0)
    SaI = inv(Sa);
else
    SaI = diag(Sa);
    SaI = diag(1./SaI);
end
SoI = diag(So);
SoI = diag(1./SoI);

%%% Invert with the n-form from Rodgers
gamma = LM_param.gamma;
%fprintf('Matrix Dimensions are as follows:')
%fprintf('size of SaI=')
%size(SaI)
%fprintf('Size of SoI=')
%size(SoI)
%fprintf('Size of k = ')
%size(K)
%fprintf('Size of gamma=')
%fprintf('Size of gamma=')
%size(gamma)
%fprintf('size of F=')
%size(F)
%fprintf('size of y=')
%size(y)
%fprintf('size of xi=')
%size(xi)
%fprintf('size of xa=')
%size(xa)
% Write out some diagnostics to see what was wrong
matr.SaI = SaI;
matr.K = K;
matr.gamma = gamma;
matr.SoI = SoI;
matr.y = y;
matr.F = F;
matr.xi = xi;
matr.xa = xa;
matrx.Sa = Sa;

LHS   = (SaI + K'*SoI*K + gamma*SaI);
RHS   = (K'*SoI * (y - F) - SaI*(xi - xa));
dx    = LHS \ RHS;
xhat  = xi + dx;

%%% Get the posterior covariance matrix and the averaging kernel matrix
S_hat  = inv(K'*SoI*K + SaI);
AK_hat = eye(size(Sa)) - S_hat*SaI;
C_hat  = sqrt(diag(S_hat));
C_hat  = diag(1./C_hat);
C_hat  = C_hat * S_hat * C_hat;

%%% Accept this step?
chi2_old = LM_param.chi2;
chi2_new = (xhat - xa)' / (SaI + K'*SoI*K) * (xhat - xa);
if ~isnan(chi2_old) && gamma > 0
    if chi2_new < chi2_old
        gamma = gamma / 3;
    else
        gamma    = gamma * 2;
        xhat     = xi;
        chi2_new = chi2_old;
    end
end

%%% Store the Levenberg-Marquardt parameters
LM_param.gamma = gamma;
LM_param.chi2  = chi2_new;

%%% Disassemble the state vector
[ems_hat, IC_hat] = disassembleStateVector(xhat,nT,nE,nI);

%%% Make the other output
out = {ems_hat,IC_hat,S_hat,AK_hat,C_hat,Sa};

end

%%% Fill diagonals
function [ sig ] = fillDiagonalsAnal(sigO,tau,St)

%%% Populate the diagonals of the covariance matrices (temporal correlations)
sig = sigO;
for k = 1:length(tau)
    if tau(k) > 0
        kk = (k-1)*length(St) + 1;
        for i = kk:(kk+length(St)-1)
            ii = mod(i-1,length(St)) + 1;
            for j = kk:(kk+length(St)-1)
                jj = mod(j-1,length(St)) + 1;
                sig(i,j) = sqrt(sigO(i,i))*sqrt(sigO(j,j)) * exp(-abs(St(ii) - St(jj))/tau(k));
            end
        end
    end
end

end


%%% =======================================================================
%%% = END
%%% =======================================================================

