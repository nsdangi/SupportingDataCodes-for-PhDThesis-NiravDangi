function [Su,Sv,Sw,Suw,Svw,Suv] = func_KaimalModel_USABL(U,Z,f,u_star,vhubplease,sigu,L_mo,ZI_in)
% Inputs:
% U: matrix [Ny x Nz] of mean wind velocity (in m/s) at each node of a grid.
% Z: matrix [Ny x Nz] of height (in m) at each node of a grid.
% f: vector [1 x Nfreq] of frequency (in Hz)
% u_star: scalar [1 x 1] friction velocity (in m/s)
% vhubplease : mean of u component at hub height (m/s)
% sigu : std of u component at hub height (m/s)
% L_MO : Monin-Obukhov length (m)
% ZI_in : Inversion layer height (m)
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
ZI = ZI_in;            % mixing layer depth (input)
L = L_mo;               % Monin–Obukhov length (positive for stable)

for jj = 1:Nm
    % --- Local variables ---
    ubar =(dummyU(jj));     % mean velocity at height z
    ubarOG=ubar;
    if ubar<=0
        ubar=abs(ubar);
    end
    zloc = dummyZ(jj);     % local height
    Ustar = u_star;        % friction velocity (scalar or array)
   
   
    
    % --- Su (streamwise component) ---
    term1_u = 0.5 .* (ZI ./ ubar) .* (ZI ./ (-L)).^(2/3);
    denom1_u = (1 + 2.2 .* (f .* ZI ./ ubar).^(5/3));
    
    term2_u = 105 .* (zloc ./ (ubar .* (1 + 15 .* zloc ./ ZI))) .* (1 - zloc ./ ZI).^2;
    denom2_u = (1 + 33 .* f .* zloc ./ (ubar .* (1 + 15 .* zloc ./ ZI))).^(5/3) .* ...
                (1 + 15 .* zloc ./ ZI).^(2/3);
    
    Su(jj,:) = Ustar.^2 .* ( term1_u ./ denom1_u + term2_u ./ denom2_u );
    
    
    % --- Sv (lateral component) ---
    term1_v = 0.95 .* (ZI ./ ubar) .* (ZI ./ (-L)).^(2/3);
    denom1_v = (1 + 2 .* f .* ZI ./ ubar).^(5/3);
    
    term2_v = 17 .* (zloc ./ (ubar .* (1 + 2.8 .* zloc ./ ZI))) .* (1 - zloc ./ ZI).^2;
    denom2_v = (1 + 9.5 .* f .* zloc ./ (ubar .* (1 + 2.8 .* zloc ./ ZI))).^(5/3) .* ...
                (1 + 2.8 .* zloc ./ ZI).^(2/3);
    
    Sv(jj,:) = Ustar.^2 .* ( term1_v ./ denom1_v + term2_v ./ denom2_v );
    
    
    % --- Sw (vertical component) ---
    term1_w = 0.95 .* (ZI ./ ubar) .* (ZI ./ (-L)).^(2/3);
    denom1_w = (1 + 2 .* f .* ZI ./ ubar).^(5/3);
    
    shape_w = (( (f .* zloc ./ ubar).^2 + (0.3 .* zloc ./ ZI).^2 ) ./ ...
              ( (f .* zloc ./ ubar).^2 + (0.15).^2 )).^(1/2);
    
    term2_w = 2 .* (zloc ./ ubar) .* (1 - zloc ./ ZI).^2;
    denom2_w = (1 + 5.3 .* (f .* zloc ./ ubar).^(5/3));
    
    Sw(jj,:) = Ustar.^2 .* ( term1_w ./ denom1_w .* shape_w + term2_w ./ denom2_w );

    if ubarOG<=0
        Su(jj,:)=-1.*Su(jj,:);
        Sv(jj,:)=-1.*Sv(jj,:);
        Sw(jj,:)=-1.*Sw(jj,:);
    end
end





end
