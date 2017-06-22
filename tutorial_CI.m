% Purpose: This script runs through how to calculate confidence intervals
%          for SHOW model estimates.
%
% Functions of interest:
%   1. FIT_SHOW_NLS.m
%   2. FIT_SHOW_LP.m
%   3. FIT_SHOW_MLE.m
%   4. OuterProduct.m
%   5. PSD_fisher_obs_LP.m
%   6. PSD_fisher_obs_MLE.m
%   7. kgradSHOW.m
%   8. PSD_gradVec_NLS.m
%   9. PSD_hessVec_NLS.m

addpath('./functions')
addpath('./data')

%% Generate PSD Time Series
rng(10)                 % Fix seed for repeatability
T_s = 5;                % Total time
SF_s = 1e6;             % Sampling frequency
f0_s = 3.5e4;           % Hz
Q_s  = 100;             % Quality factor
k_s  = 0.172;           % N/m
Kb = 1.381e-23;         % Boltzmann's constant
T = 298;                % Kelvin
CONST = 1e30;           % Unit conversion
As_s = 4*Kb*T/(k_s*Q_s*f0_s*2*pi) * CONST;  % SHO
Aw_s = 25000;                               % White noise

% Transfer Functions
f  = linspace(1/T_s,SF_s,SF_s*T_s); % freq space.
abs_C = abs((1/k_s) * 1./ ( (1-(f/f0_s).^2) + 1j*f/f0_s/Q_s));
g = k_s/2/pi/f0_s/Q_s;

xAPSD = f; % Change naming convention
yAPSD = (4*Kb*T*g) .* abs_C.^2 * CONST; % Thermal physics
yAPSD = yAPSD + Aw_s;

% TIME DOMAIN
[xTime, yTime] = tsSim(SF_s, SF_s*T_s, xAPSD, yAPSD);

% fft(yTime)
N = length(yTime)-rem(length(yTime),2); % use even PSD's only
FR = SF_s/N;
yTime = yTime(1:N);
yFreq = fft(yTime)/(sqrt(FR)*N);
% Remove DC component, select first frequency up to the Nyquist frequency.
% Drop negative frequencies. For a real signal, negative frequencies are redundant.
yFreq = yFreq(2:N/2+1);
xFreq = linspace(1/T_s, SF_s/2, length(yFreq));

yFreq(end) = real(yFreq(end)); % Nyquist freq should be real for real signal
yFreq = abs(yFreq).^2; % Take the squared magnitude at each frequency
% Make the PSD single-sided by doubling everything except the Nyquist frequency.
% This corrects for the loss in amplitude caused by dropping negative frequencies
% the Nyquist frequency at N/2 doesn't undergo folding, and therefore is NOT doubled.
yFreq(1:N/2-1) = 2.*yFreq(1:N/2-1);


%% Fit NLS, LP, and MLE models to data
PSD_x = xFreq;
PSD_y = yFreq;
B = 100; % Bin size for NLS and LP
[As_NLS, Aw_NLS, f0_NLS, Q_NLS, exitflag_NLS] = FIT_SHOW_NLS(PSD_x, PSD_y, Q_s, f0_s, Aw_s, As_s, B);
[As_LP, Aw_LP, f0_LP, Q_LP, exitflag_LP] = FIT_SHOW_LP(PSD_x, PSD_y, Q_s, f0_s, Aw_s, As_s, B);
[As_MLE, Aw_MLE, f0_MLE, Q_MLE, exitflag_MLE] = FIT_SHOW_MLE(PSD_x, PSD_y, Q_s, f0_s, Aw_s, As_s, 1);

% Back-out k estimates
k_NLS = 4*Kb*T/(As_NLS*Q_NLS*f0_NLS*2*pi) * CONST;
k_LP = 4*Kb*T/(As_LP*Q_LP*f0_LP*2*pi) * CONST;
k_MLE = 4*Kb*T/(As_MLE*Q_MLE*f0_MLE*2*pi) * CONST;

%% Fit results
T_fit = table(categorical({'NLS';'LP';'MLE'}),[As_NLS;As_LP;As_MLE],...
    [Q_NLS;Q_LP;Q_MLE],[f0_NLS;f0_LP;f0_MLE],[Aw_NLS;Aw_LP;Aw_MLE],[k_NLS;k_LP;k_MLE],...
    'VariableNames',{'Method','As' 'Q' 'f0' 'Aw' 'k'})

%% Confidence interval calculations
% Bandpass filter on frequency range (usually denoising frequency range)
freqRange_fit = [f0_s-f0_s/sqrt(2) f0_s+f0_s/sqrt(2)];
cond = PSD_x > freqRange_fit(1) & PSD_x < freqRange_fit(2); % filtering condition
xFreq = PSD_x(cond);
yFreq = PSD_y(cond);

% Data Compression (Mean-Decimation)
N = length(xFreq);
N = N-rem(N,B);
xFreq_bin = mean(reshape(xFreq(1:N), B, N/B),1); % Take average freq of each bin
yFreq_bin = mean(reshape(yFreq(1:N), B, N/B),1);  % Mean decimation

% Parameter estimates
theta_mle = [T_fit.f0(3) T_fit.Q(3) T_fit.k(3) T_fit.Aw(3)]; % parameter estimates [f0 Q k Aw]
theta_lp = [T_fit.f0(2) T_fit.Q(2) T_fit.k(2) T_fit.Aw(2)]; % parameter estimates [f0 Q k Aw]
theta_nls = [T_fit.f0(1) T_fit.Q(1) T_fit.k(1) T_fit.Aw(1)]; % parameter estimates [f0 Q k Aw]

% NLS, LP, and MLE likelihoods
bias_const = psi(B) - log(B); % needed for correction to LP
S_func = @(x) SHOW_PSD(xFreq_bin, 4*Kb*T/(x(3)*x(2)*x(1)*2*pi)*CONST, x(4), x(1), x(2));
true_nls_LL = @(x) sum(-(yFreq_bin - S_func(x)).^2./2);
true_lp_LL = @(x) sum((B/2) .* (log(yFreq_bin) - bias_const - log(SHOW_PSD(xFreq_bin, 4*Kb*T/(x(3)*x(2)*x(1)*2*pi)*CONST, x(4), x(1), x(2)))).^2);
true_mle_LL = @(x) -sum(-yFreq./SHOW_PSD(xFreq, 4*Kb*T/(x(3)*x(2)*x(1)*2*pi)*CONST, x(4), x(1), x(2)) - log(SHOW_PSD(xFreq, 4*Kb*T/(x(3)*x(2)*x(1)*2*pi)*CONST, x(4), x(1), x(2))));

% MLE hessian
gradS = kgradSHOW(xFreq,theta_mle(1),theta_mle(2),theta_mle(3));
S = SHOW_PSD(xFreq, 4*Kb*T/(theta_mle(3)*theta_mle(2)*theta_mle(1)*2*pi)*CONST, theta_mle(4), theta_mle(1), theta_mle(2));
hessS = khessSHOW(xFreq,theta_mle(1),theta_mle(2),theta_mle(3));
mle_hess = PSD_fisher_obs_MLE(yFreq', S', gradS, hessS);

% LP hessian
gradS = kgradSHOW(xFreq_bin,theta_lp(1),theta_lp(2),theta_lp(3));
S = SHOW_PSD(xFreq_bin, 4*Kb*T/(theta_lp(3)*theta_lp(2)*theta_lp(1)*2*pi)*CONST, theta_lp(4), theta_lp(1), theta_lp(2));
z = log(yFreq_bin) - bias_const;
hessS = khessSHOW(xFreq_bin,theta_lp(1),theta_lp(2),theta_lp(3));
lp_hess = PSD_fisher_obs_LP(S, gradS, hessS, z, B);

% NLS sandwich estimator
S = SHOW_PSD(xFreq_bin, 4*Kb*T/(theta_nls(3)*theta_nls(2)*theta_nls(1)*2*pi)*CONST, theta_nls(4), theta_nls(1), theta_nls(2));
gradS = kgradSHOW(xFreq_bin,theta_nls(1),theta_nls(2),theta_nls(3));
hessS = khessSHOW(xFreq_bin,theta_nls(1),theta_nls(2),theta_nls(3));
[B,myGrad] = PSD_gradVec_NLS(S, gradS, yFreq_bin);
A = PSD_hessVec_NLS(S, gradS, hessS, yFreq_bin);
V = inv(-A)' * B * inv(-A); % Sandwich estimator matrix

% 95% confidence wintervals, parameter estimates [f0 Q k Aw]
se_mle = sqrt(diag(inv(mle_hess)))';
se_lp = sqrt(diag(inv(lp_hess)))';
se_nls = sqrt(diag(V))';

nlsCI = table(categorical({'f0';'Q';'k'}),...
    [T_fit.f0(1);T_fit.Q(1);T_fit.k(1)],...
    [se_nls(1);se_nls(2);se_nls(3)],...
    'VariableNames',{'Parameter','Estimate','SE'})

lpCI = table(categorical({'f0';'Q';'k'}),...
    [T_fit.f0(2);T_fit.Q(2);T_fit.k(2)],...
    [se_lp(1);se_lp(2);se_lp(3)],...
    'VariableNames',{'Parameter','Estimate','SE'})

mleCI = table(categorical({'f0';'Q';'k'}),...
    [T_fit.f0(3);T_fit.Q(3);T_fit.k(3)],...
    [se_mle(1);se_mle(2);se_mle(3)],...
    'VariableNames',{'Parameter','Estimate','SE'})
