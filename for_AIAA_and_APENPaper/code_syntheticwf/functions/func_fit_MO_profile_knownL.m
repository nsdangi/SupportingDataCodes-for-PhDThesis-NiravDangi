function params = func_fit_MO_profile_knownL(z, u_obs, z0, L)
% FIT_MO_PROFILE_KNOWNL
% Fit only ustar, L is fixed

    % Initial guess
    ustar0 = 0.3;

    % Cost function
    fun = @(ustar) mean((func_wind_profile_MO(z, z0, ustar, L) - u_obs).^2);

    % Optimization
    options = optimset('Display','off');
    ustar_opt = fminsearch(fun, ustar0, options);

    params.ustar = ustar_opt;
    params.L     = L;   % just return input
end