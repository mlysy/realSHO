% Purpose: This script runs through how to fit for SHOW estimates, baseline
%          enviroment WITH added sine-wave noise.
%
% Data of interest:
%   1. ./Data/Q_100_PSD_noise.mat
%
% Functions of interest:
%   1. FIT_SHOW_NLS.m
%   2. FIT_SHOW_LP.m
%   3. FIT_SHOW_MLE.m
%   4. psd_denoise.m
%   5. fisherGstat.m

addpath('./functions')

%% Load periodgram data
load('./data/Q_100_PSD_noise.mat')
PSD_x = Q_100_PSD_noise.xPSD;
PSD_y = Q_100_PSD_noise.yPSD;

%% True Simulation Parameters
B = 100;                % Size of bins for fitting
f0_s = 33553;           % Resonance Frequency (Hz)
Q_s  = 100;             % Quality factor
k_s  = 0.172;           % Stiffness (N/m)
Kb = 1.381e-23;         % Boltzmann's constant
T = 298;                % Kelvin
CONST = 1e30;           % Unit conversion
As_s = 4*Kb*T/(k_s*Q_s*f0_s*2*pi) * CONST;  % SHO "amplitude"
Aw_s = 19000;                               % white noise amplitude

%% Fit NLS, LP, and MLE models to data
% First, remove systematic noise from PSD
pcut = 0.01; % p-value to use in periodicity removal
freq_range = [f0_s-f0_s/sqrt(2) f0_s+f0_s/sqrt(2)]; % Frequency range to denoise
[clean_xPSD, clean_yPSD] = psd_denoise(PSD_x, PSD_y, Q_s, f0_s, Aw_s, As_s, freq_range, pcut);

% Now fit
[As_NLS, Aw_NLS, f0_NLS, Q_NLS, exitflag_NLS] = FIT_SHOW_NLS(clean_xPSD, clean_yPSD, Q_s, f0_s, Aw_s, As_s, B);
[As_LP, Aw_LP, f0_LP, Q_LP, exitflag_LP] = FIT_SHOW_LP(clean_xPSD, clean_yPSD, Q_s, f0_s, Aw_s, As_s, B);
[As_MLE, Aw_MLE, f0_MLE, Q_MLE, exitflag_MLE] = FIT_SHOW_MLE(clean_xPSD, clean_yPSD, Q_s, f0_s, Aw_s, As_s, 1);

% Back-out k estimates
k_NLS = 4*Kb*T/(As_NLS*Q_NLS*f0_NLS*2*pi) * CONST;
k_LP = 4*Kb*T/(As_LP*Q_LP*f0_LP*2*pi) * CONST;
k_MLE = 4*Kb*T/(As_MLE*Q_MLE*f0_MLE*2*pi) * CONST;

%% Fit results
T_fit = table(categorical({'NLS';'LP';'MLE'}),[As_NLS;As_LP;As_MLE],...
    [Q_NLS;Q_LP;Q_MLE],[f0_NLS;f0_LP;f0_MLE],[Aw_NLS;Aw_LP;Aw_MLE],[k_NLS;k_LP;k_MLE],...
    'VariableNames',{'Method','As' 'Q' 'f0' 'Aw' 'k'})
