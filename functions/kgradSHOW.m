function [grad] = kgradSHOW(f,f0,Q,k)
% Creation Date: Nov. 27, 2015
% Purpose: Calculates gradient of SHO PSD model from AFM data
% Model: SHOW = Aw + As/((1-(f/f0)^2)^2 + (f/Q/f0)^2)
% Input:
%   f:  frequency
%   f0: resonance frequency
%   Q:  quality factor
%   
% Output: Matrix [dSHO/df0 dSHO/dQ dSHO/dk]

% Physcial Constants
Kb = 1.381e-23; % Boltzmann's constant
T = 298;        % Kelvin
CONST = 1e30;   % yes, unit conversion

rr = (f./f0).^2 - 1;
rr2 = ((f./f0).^2 - 1).^2;

PHI = k.*Q.*f0.*rr2 + (f.^2).*k./f0./Q;
PSI = -2.*Kb.*T.*CONST./pi./(PHI.^2);

dPHI_df0 = k.*Q.*rr2 - 4.*k.*Q.*f0.*rr.*f.^2./(f0^3) - (f./f0).^2.*(k/Q);
dPHI_dQ = k.*f0.*rr2 - (f.^2)./f0 .* k./(Q^2);
dPHI_dk = Q.*f0.*rr2 + (f.^2)./f0./Q;

dS_df0 = PSI .* dPHI_df0;
dS_dQ = PSI .* dPHI_dQ;
dS_dk = PSI .* dPHI_dk;
dS_dAw = ones(1,length(f));

% Ensure correct output (convert all to column vectors)
dS_df0 = reshape(dS_df0, [max(size(dS_df0)) 1]);
dS_dQ = reshape(dS_dQ, [max(size(dS_dQ)) 1]);
dS_dk = reshape(dS_dk, [max(size(dS_dk)) 1]);
dS_dAw = reshape(dS_dAw, [max(size(dS_dAw)) 1]);

% Gradient
grad = [dS_df0 dS_dQ dS_dk dS_dAw];
end

