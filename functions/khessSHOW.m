function [hess] = khessSHOW(f,f0,Q,k)
% Creation Date: Nov. 27, 2015
% Purpose: Calculates Hessian of SHO PSD model from AFM data
% Model: SHOW = Aw + As/((1-(f/f0)^2)^2 + (f/Q/f0)^2)
% Input:
%   f:  frequency
%   f0: resonance frequency
%   Q:  quality factor
%   k:  Stiffness
%   Aw: Amplitude
% Output: M x 4 x 4 where M is length of input frequency vector f,
% so output "slice" is of form:
% [d2SHO/df02   d2SHO/df0dQ d2SHO/df0dk d2SHO/df0dAw
%  d2SHO/dQdf0  d2SHO/dQ2   d2SHO/dQdk d2SHO/dQdAw
%  d2SHO/dAsdf0 d2SHO/dAsdQ d2SHO/dk2 d2SHO/dAsdAw
%  d2SHO/dAwdf0 d2SHO/dAwdQ d2SHO/dAwdk d2SHO/dAw2]

% Physcial Constants
Kb = 1.381e-23; % Boltzmann's constant
T = 298;        % Kelvin
CONST = 1e30;   % yes, unit conversion

rr = (f./f0).^2 - 1; % Good
rr2 = ((f./f0).^2 - 1).^2; % Good

PHI = k.*Q.*f0.*rr2 + (f.^2).*k./f0./Q; % Good
PSI = -2.*Kb.*T.*CONST./pi./(PHI.^2); % Good

% Derivatives of PHI
dPHI_df0 = k.*Q.*rr2 - 4.*k.*Q.*f0.*rr.*f.^2./(f0^3) - (f./f0).^2.*(k/Q); % Good
dPHI_dQ = k.*f0.*rr2 - (f.^2)./f0 .* k./(Q^2); % Good
dPHI_dk = Q.*f0.*rr2 + (f.^2)./f0./Q; % Good

dPHI2_df02 = -4.*k.*Q.*rr.*f.^2./f0^3 + ...
            8.*k.*Q.*f.^2./f0^3.*rr2 + ...
            16.*k.*Q.*(f./f0).^2.*rr.*f.^2./f0^3 - ...
            2.*f.^2./f0^3.*k./Q; % Done
dPHI2_dQ2 = 2.*f.^2./f0.*k/Q.^3; % Done
dPHI2_dk2 = 0; % Done
dPHI2_dQdf0 = k.*rr2 - 4.*k.*f0.*rr.*f.^2./(f0^3) + (f./f0).^2.*(k/Q^2); % Done
dPHI2_dkdf0 = Q.*rr2 - 4.*Q.*f0.*rr.*f.^2./(f0^3) - (f./f0).^2./Q; % Done
dPHI2_dkdQ = f0.*rr2 + (f.^2)./f0./(Q^2); % Done

% Derivatives of PSI
dPSI_df0 = 4.*Kb.*T.*CONST./pi./(PHI.^3).*dPHI_df0; % Good
dPSI_dQ = 4.*Kb.*T.*CONST./pi./(PHI.^3).*dPHI_dQ; % Good
dPSI_dk = 4.*Kb.*T.*CONST./pi./(PHI.^3).*dPHI_dk; % Good

% Hessian components
dS2_df02 = dPSI_df0 .* dPHI_df0 + PSI .* dPHI2_df02; % Good
dS2_dQ2 = dPSI_dQ .* dPHI_dQ + PSI .* dPHI2_dQ2; % Good
dS2_dk2 = dPSI_dk .* dPHI_dk + PSI .* dPHI2_dk2; % Good
dS2_dkdf0 = dPSI_dk .* dPHI_df0 + PSI .* dPHI2_dkdf0; % Good
dS2_dQdf0 = dPSI_dQ .* dPHI_df0 + PSI .* dPHI2_dQdf0; % Good
dS2_dkdQ = dPSI_dk .* dPHI_dQ + PSI .* dPHI2_dkdQ; % Good

% Hessian M x 3 x 3
M = length(f);
hess = zeros(M,4,4);
hess(:,1,1) = dS2_df02;
hess(:,1,2) = dS2_dQdf0;
hess(:,1,3) = dS2_dkdf0;
hess(:,2,1) = dS2_dQdf0;
hess(:,2,2) = dS2_dQ2;
hess(:,2,3) = dS2_dkdQ;
hess(:,3,1) = dS2_dkdf0;
hess(:,3,2) = dS2_dkdQ;
hess(:,3,3) = dS2_dk2;

hess(:,1,4) = 0;
hess(:,2,4) = 0;
hess(:,3,4) = 0;
hess(:,4,4) = 0;
hess(:,4,1) = 0;
hess(:,4,2) = 0;
hess(:,4,3) = 0;
hess(:,4,4) = 0;
end

