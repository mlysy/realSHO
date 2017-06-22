function [As, Aw, f0, Q, exitflag] = FIT_SHOW_MLE(PSD_x, PSD_y, Q_s, f0_s, Aw_s, As_s, bin_size, varargin)
% Performs Whittle log-likelihood based maximum likelihood estimation for
% SHOW (simple harmonic oscillator + white noise) model. Profile likelihood
% is used.
%
% Inputs:
%   PSD_x: Periodogram frequency values (no binning)
%   PSD_y: Periodogram ordinates (no binning)
%     Q_s: Initial quality factor estimate
%    f0_s: Initial resonance frequency estimate
%    Aw_s: Initial white noise amplitude estimate
%    As_s: Initial SHO "amplitude" estimate
%    bin_size: chosen bin size to use for fitting (integer>=1)
%    Varargin: If last argument is missing, fitting frequency range is auto
%              set to f0_s +/- f0_s/sqrt(2)
%
% Outputs:
%    As: SHO "amplitude" estimate
%    Aw: White noise amplitude estimate
%    f0: Resonance frequency estimate
%     Q: Quality factor estimate
%  exitflag: 0 = "success", 1 = "fminsearch() convergence issues"
%
% Example Usage:
% [As, Aw, f0, Q, exitflag] = FIT_SHOW_MLE(PSD_x, PSD_y, Q_s, f0_s, Aw_s, As_s, bin_size)

% Output parameter init
As = NaN;
Aw = NaN; 
f0 = NaN;
Q = NaN;
exitflag = 1; % 1 = fminsearch failed, 0 = success

%% --- DETERMINE RANGE OF DATA TO FIT -------------------------------------
if length(varargin) < 1
    f_lowerBound = f0_s - f0_s/sqrt(2);
    f_upperBound = f0_s + f0_s/sqrt(2);
else
    f_lowerBound = varargin{1}(1);
    f_upperBound = varargin{1}(2);
end
fprintf('Frequency Cut-Off: Lower is %d, Upper is %d \n',f_lowerBound,f_upperBound)
cond = PSD_x > f_lowerBound & PSD_x < f_upperBound; % Filtering condition
xPSD = PSD_x(cond);
yPSD = PSD_y(cond);

%% Data Compression (Mean-Decimation)
f2 = xPSD;      % Vector of frequencies
N = length(f2);
N = N-rem(N,bin_size);
f2 = mean(reshape(f2(1:N), bin_size, N/bin_size),1).^2;   % Take average freq of each bin
S = mean(reshape(yPSD(1:N), bin_size, N/bin_size),1);     % Mean decimation

%% SHOW Parameter Fit
f0 = f0_s;               % Resonance frequency
Q = Q_s;                 % Quality factor
Arw = Aw_s/As_s;     	 % Relative Aw (Arw = Aw/As)

% ----First fit for Q----
[SHOFWfit,~,exitflag] = fminsearch(@(x) TrueMLEObjFun(S, f2, Arw, f0, x(1)), [Q]);
if exitflag ~= 1
    disp('fminsearch failed to converge')
    return
end
Q = SHOFWfit(1);

% ----Now do f0 and Q together----
[SHOFWfit,~,exitflag] = fminsearch(@(x) TrueMLEObjFun(S, f2, Arw, x(1), x(2)), [f0 Q]);
if exitflag ~= 1
    disp('fminsearch failed to converge')
    return
end
f0 = SHOFWfit(1);
Q = SHOFWfit(2);

% ----Now do Aw, f0, and Q together----
[SHOFWfit,~,exitflag] = fminsearch(@(x) TrueMLEObjFun(S, f2, x(1), x(2), x(3)), [Arw f0 Q]);
if exitflag ~= 1
    disp('fminsearch failed to converge')
    return
end
Arw = SHOFWfit(1);
f0 = SHOFWfit(2);
Q = SHOFWfit(3);

% Final estimates and return
As = mean(S./G_PSD(sqrt(f2), Arw, f0, Q)); % Back-out As from profile parameter
Aw = As*Arw; % Backout the Aw from Arw
end

function F = TrueMLEObjFun(y, f2, Arw, f0, Q)
S = G_PSD(sqrt(f2), Arw, f0, Q);
tau = mean(y./S);
F = -PSD_loglik_PL(y, S, tau);
end

function S = G_PSD(f, Arw, f0, Q)
S = Arw + 1./((1 - (f/f0).^2).^2 + (f/Q/f0).^2);
end

function ll = PSD_loglik_PL(y, S, tau)
ll = -sum(y./(tau.*S) + log(tau) + log(S));
end

