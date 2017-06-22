function S = SHOW_PSD(f, As, Aw, f0, Q)
% Simple Harmonic Oscillator + White Noise model. Returns PSD ordinates at
% each frequency value in f.
%
% Inputs:
%    f: Periodogram frequency values (no binning)
%   As: SHO "amplitude"
%   Aw: White noise amplitude
%   f0: Resonance frequency
%    Q: Quality factor estimate
%
% Outputs:
%    S: Power Spectral Density
%
% Example Usage:
% S = SHOW_PSD(f, As, Aw, f0, Q)

S = Aw + As .* 1./((1 - (f/f0).^2).^2 + (f/Q/f0).^2);
end

