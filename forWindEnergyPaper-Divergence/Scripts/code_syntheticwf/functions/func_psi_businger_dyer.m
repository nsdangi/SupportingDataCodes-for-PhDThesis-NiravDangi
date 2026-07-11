function psi = func_psi_businger_dyer(z, L)
% Standard BD formulation with Högström constants

    gamma = 19.3;

    x = (1 - gamma * (z ./ L)).^(1/4);

    psi = 2 * log((1 + x) / 2) ...
        + log((1 + x.^2) / 2) ...
        - 2 * atan(x) ...
        + pi/2;
end