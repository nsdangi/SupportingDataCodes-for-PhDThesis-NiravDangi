function u = func_wind_profile_MO(z, z0, ustar, L, kappa)
% WIND_PROFILE_MO
% Compute Monin-Obukhov wind profile using provided L
%   Stable   -> Holtslag
%   Unstable -> Businger-Dyer
%   Neutral  -> Log-law
    if nargin < 5
        kappa = 0.4;
    end

    % Neutral case
    if abs(L) > 1e5
        psi_z  = 0;
        psi_z0 = 0;
    else
        psi_z  = func_psi_combined(z, L);
        psi_z0 = func_psi_combined(z0, L);
    end

    u = (ustar / kappa) .* (log(z ./ z0) - psi_z + psi_z0);
end