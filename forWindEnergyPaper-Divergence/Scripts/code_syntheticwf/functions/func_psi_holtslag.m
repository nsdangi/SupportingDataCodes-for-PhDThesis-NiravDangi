function psi = func_psi_holtslag(z, L)

    zL = z ./ L;

    psi = - zL ...
          - (2/3) * (zL - (5/0.35)) .* exp(-0.35 * zL) ...
          - (10 / 1.05);
end