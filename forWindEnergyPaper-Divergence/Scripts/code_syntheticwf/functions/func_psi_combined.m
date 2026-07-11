function psi = func_psi_combined(z, L)
% PSI_COMBINED
% Switch between:
%   Stable   -> Holtslag
%   Unstable -> Businger-Dyer

    if L > 0
        psi = func_psi_holtslag(z, L);
    elseif L < 0
        psi = func_psi_businger_dyer(z, L);
    else
        psi = 0;
    end
end