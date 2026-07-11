function [Su,Sv,Sw,Suw,Svw,Suv] = func_KaimalModel_IECKAI(U,Z,f,u_star,vhubplease,sigu)
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
for jj=1:Nm
    fr = (f*dummyZ(jj)/dummyU(jj));
    %    Su(jj,:) = 102.*fr./(1+33.*fr).^(5/3).*u_star.^2./f; % Kaimal  model (NOT normalized)
    Su(jj,:) = (4.*f.*340.2)./(vhubplease.*(1+(6.*f.*340.2)./(vhubplease)).^(5/3)).*(((sigu).^2)./(f));
    %    Sv(jj,:) = 17.*fr./(1+9.5.*fr).^(5/3).*u_star.^2./f; % Kaimal  model (NOT normalized)
    Sv(jj,:) = (4.*f.*113.4)./(vhubplease.*(1+(6.*f.*113.4)./(vhubplease)).^(5/3)).*(((0.8.*sigu).^2)./(f));
    %    Sw(jj,:) = (2.*fr./(1+5.*fr.^(5/3))).*u_star.^2./f; % Kaimal model (NOT normalized)
    Sw(jj,:) = (4.*f.*27.72)./(vhubplease.*(1+(6.*f.*27.72)./(vhubplease)).^(5/3)).*(((0.5.*sigu).^2)./(f));
    Suw(jj,:) = -14.*fr./(1+10.5*fr).^(7/3).*(0.5*u_star).^2./f; % Kaimal cross-spectrum model (NOT normalized)
    Svw(jj,:) = -14.*fr./(1+10.5*fr).^(7/3).*(0.5*u_star).^2./f; % Kaimal cross-spectrum model (NOT normalized)
    Suv(jj,:) = 14.*fr./(1+10.5*fr).^(7/3).*(0.25*u_star).^2./f; % Kaimal cross-spectrum model (NOT normalized)
end





end