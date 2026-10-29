function [Su,Sv,Sw,Suw,Svw,Suv] = func_KaimalModel_NABL(U,Z,f,u_star,vhubplease,sigu)
% Inputs:
% U: matrix [Ny x Nz] of mean wind velocity (in m/s) at each node of a grid.
% Z: matrix [Ny x Nz] of height (in m) at each node of a grid.
% f: vector [1 x Nfreq] of frequency (in Hz)
% u_star: scalar [1 x 1] friction velocity (in m/s)
% vhubplease : mean of u component at hub height (m/s)
% sigu : std of u component at hub height (m/s)
% Outputs:
% Posing Nm = Nx*Ny and recalling that PSD = power spectral density and
% CPSD = cross-power spectral density:
% Su: vector [Nm x 1] corresponding to the PSD the u-component
% Sv: vector [Nm x 1] corresponding to the PSD the v-component.
% Sw: vector [Nm x 1] corresponding to the PSD the w-component
% Suw: vector [Nm x 1] corresponding to the CPSD the u and w components
%
% Original Author: E. Cheynet - uiB - Last modified 2024-04-16

% ##### Modified version ##### %

%%
Nfreq = numel(f);
Nm = numel(U(:));
Su = zeros(Nm,Nfreq); % preallocation
Sv = zeros(Nm,Nfreq); % preallocation
Sw = zeros(Nm,Nfreq); % preallocation
Suw = zeros(Nm,Nfreq); % preallocation
Svw = zeros(Nm,Nfreq); % preallocation
Suv = zeros(Nm,Nfreq); % preallocation
dummyU = U(:);
dummyZ = Z(:);
L_mo=10^10;
L = L_mo;               % Monin–Obukhov length (positive for stable)
for jj = 1:Nm
    % --- Local variables ---
    ubar = dummyU(jj);     % mean velocity at height z
    zloc = dummyZ(jj);     % local height
    Ustar = u_star;         % friction velocity
   
    
    % --- Adjust L to keep 0 < z/L < 2 ---
    z_over_L = zloc / L;
    if z_over_L <= 0
        L = zloc / 0.01;   % small positive value to avoid division by zero
        z_over_L = 0.01;
    elseif z_over_L > 2
        L = zloc / 2;
        z_over_L = 2;
    end
    
    % --- Stability functions ---
    phiM = 1 + 4.7 * z_over_L;
    phiE = (1 + 2.5 * z_over_L^0.6)^(3/2);  % phi_E^(2/3) given
    
    % --- Dimensionless frequency ---
    fr = f .* zloc ./ (ubar .* phiM);
    
    % --- SMOOTH model coefficients ---
    s1u = 79.0;  s2u = 263.0;
    s1v = 13.0;  s2v = 32.0;
    s1w = 3.5;   s2w = 8.6;
    
    % --- Compute spectra ---
    Su(jj,:) = Ustar.^2 .* ...
        ( s1u .* (zloc ./ (ubar .* phiM)) .* ((phiE ./ phiM).^(2/3)) ) ./ ...
        ( 1.0 + s2u .* fr.^(5/3) );
    
    Sv(jj,:) = Ustar.^2 .* ...
        ( s1v .* (zloc ./ (ubar .* phiM)) .* ((phiE ./ phiM).^(2/3)) ) ./ ...
        ( 1.0 + s2v .* fr.^(5/3) );
    
    Sw(jj,:) = Ustar.^2 .* ...
        ( s1w .* (zloc ./ (ubar .* phiM)) .* ((phiE ./ phiM).^(2/3)) ) ./ ...
        ( 1.0 + s2w .* fr.^(5/3) );
end

end



