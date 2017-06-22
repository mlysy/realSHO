function [clean_xPSD, clean_yPSD] = psd_denoise(xPSD, yPSD, Q, f0, Aw, As, freq_range, pcut)
% Performs systematic noise removal using a long-standing test for spectral
% periodicites. The fisher g-statistics is calculated, and ordinates are
% replaced sequentially based on p-value cut-off pcut.
%
% Inputs:
%   xPSD: Periodogram frequency values (no binning)
%   yPSD: Periodogram ordinates (no binning)
%     Q: Initial quality factor estimate
%    f0: Initial resonance frequency estimate
%    Aw: Initial white noise amplitude estimate
%    As: Initial SHO "amplitude" estimate
%    freq_range: Frequency range to denoise
% Output:
%   clean_xPSD: Cleaned periodogram frequencies
%   clean_xPSD: Cleaned periodogram ordinates
% Example Usage:
% [clean_xPSD, clean_yPSD] = psd_denoise(xPSD, yPSD, Q, f0, Aw, As, freq_range, pcut)

%% --- DETERMINE RANGE OF DATA TO USE -------------------------------------
f_lowerBound = freq_range(1);
f_upperBound = freq_range(2);
fprintf('Frequency Cut-Off: Lower is %d, Upper is %d \n',f_lowerBound,f_upperBound)
cond = xPSD > f_lowerBound & xPSD < f_upperBound; % Filtering condition
PSD_x = xPSD(cond);
PSD_y = yPSD(cond);

%% --- REMOVE PERIODIC COMPONENTS FROM PSD --------------------------------
Arw = Aw/As;
q = length(PSD_y);
yPSD_Fit = exp(logSHOW(PSD_x, As, Arw, f0, Q));
gPSD = PSD_y./yPSD_Fit;
[a, inda] = max(gPSD/sum(gPSD));
g = fisherGstat(a,q);
ncut = 0;
while g < pcut
    fprintf('G-stat found periodicity, replacing freq %d \n', PSD_x(inda));
    ncut = ncut+1;
    PSD_y(inda) = exprnd(exp(logSHOW(PSD_x(inda), As, Arw, f0, Q)));
    gPSD = PSD_y./yPSD_Fit;
    [a, inda] = max(gPSD/sum(gPSD));
    g = fisherGstat(a,q);
end
yPSD(cond) = PSD_y; % replace cleaned component of PSD
clean_yPSD = yPSD;
clean_xPSD = xPSD;
end

function logS = logSHOW(f, As, Arw, f0, Q)
logS = log(As) + log(Arw + 1./((1 - (f/f0).^2).^2 + (f/f0/Q).^2));
end

