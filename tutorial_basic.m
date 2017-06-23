% Purpose: This script runs through how to fit for SHOW estimates, baseline
%          enviroment WITHOUT added sine-wave noise, in addition to time
%          series to PSD conversion.
%
%
% Functions of interest:
%   1. FIT_SHOW_NLS.m
%   2. FIT_SHOW_LP.m
%   3. FIT_SHOW_MLE.m
%   4. tsSim.m

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
[xFreq,yFreq] = get_periodogram(yTime, SF_s, T_s);

%% Visual in both time and frequency domain
% Perform binning top help visuals
bin_size = 100;
f2 = xFreq;
N = length(f2);
N = N-rem(N,bin_size);
f2 = mean(reshape(f2(1:N), bin_size, N/bin_size),1);   % Take average freq of each bin
S = mean(reshape(yFreq(1:N), bin_size, N/bin_size),1);  % Mean decimation

subplot(2,1,1)
plot(xTime,yTime)
legend('Time Domain')
xlabel('Time (s)')
ylabel('Signal')

subplot(2,1,2)
loglog(f2,S,xAPSD,yAPSD)
legend('Frequency Domain','Theoretical PSD','location','northwest')
xlabel('Frequency (Hz)')
ylabel('PSD')

%% SHO Fit using NLS, LP, and MLE
PSD_x = xFreq;
PSD_y = yFreq;
B = 100; % Bin size
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



