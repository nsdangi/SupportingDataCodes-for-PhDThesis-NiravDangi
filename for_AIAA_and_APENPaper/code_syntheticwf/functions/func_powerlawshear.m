function u_z = func_powerlawshear(u_ref, z_ref, z_target, alpha)
% powerLawShear Calculates wind speed at a target height using the Power Law.
%
% Inputs:
%   u_ref    - Wind speed at reference height (m/s)
%   z_ref    - Reference height (m)
%   z_target - Height(s) at which to calculate wind speed (m)
%   alpha    - Wind shear exponent (default 0.143)

u_z = u_ref .* (z_target ./ z_ref).^alpha;
end