function [xFreq,yFreq] = get_periodogram(yTime,SF_s,T_s)
% Calculates the one-sided periodogram from given time series and sampling
% frequency SF
%
% Inputs:
%   yTime: Time series to transform
%   SF_s: Sampling frequecy (Hz)
%    T_s: Total time (seconds)
%
% Outputs:
%    xFreq: One-sided periodogram frequencies 
%    yFreq: One-sidede periodogram ordinates
%
% Example Usage:
% [xFreq,yFreq] = periodogram(yTime,SF)

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
end

