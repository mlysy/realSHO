function [xTime, yTime] = tsSim(SF, N, xPSD, yPSD)
% Simulates a time series from a given PSD at a given sampling frequency
%
% Inputs:
%    SF: Sampling frequency
%     N: Total number of samples
%  xPSD: Frequency basis, possibly including DC offset (xPSD(1) = 0)
%        if DC offset is not provided then it is set to 0.
%  yPSD: Power spectrum up to Nyquist frequency
%
% Outputs:
%  xTime: Time observation times
%  yTime: Time series (i.e., amplitude)
%
% Example Usage:
%   [xTime, yTime] = tsSim(SF, N, xPSD, yPSD)
  
% Sampling frequency
N = 2*N; % Double time series length, throw out half the data later

% Error handling
if xPSD(end)*2 < SF
    error('Sampling rate exceeds the range of the spectrum.');
end
if 1/xPSD(1+ (xPSD(1) == 0)) < (N-1)/SF
    disp(sprintf('%s\n%s', 'Length of time series exceeds spectral period.', ...
        'Lower frequencies taken to be minimum frequency value.'));
end

m = floor(N/2);
fTs = (1:m)*(SF/N); % Frequency range of the time series
PTs = interp1(xPSD, yPSD, fTs, 'linear', yPSD(1)); % Spectrum of the time series

% Simulate a Fourier basis
yFreq = randn(size(fTs)) + 1i*randn(size(fTs));

% One sqrt(2) is for unfolding, the other is for the mean of the sum of
% squares of two iid Gaussians
yFreq = yFreq.* sqrt(PTs)/2;

% Unfold (Nyquist frequency doesn't undergo aliasing and is real)
yFreq(end) = real(yFreq(end))*sqrt(2);
yConj = conj(wrev(yFreq(1:end-1)));

dcOffset = (xPSD(1) == 0)*yPSD(1); % DC offset

yFreq = [dcOffset yFreq yConj];

% Recover the time series
FR = SF/N;
yTime = ifft(yFreq*sqrt(FR)*N);
yTime = real(yTime(1:end/2));
xTime = (0:length(yTime)-1)/SF;
end

