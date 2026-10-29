function [u_corr,v_corr,w_corr] = func_windSimFaster_NABL(Y,Z,U,V,W,f,u,v,w)
% windSimFaster Introduces spatial correlation into uncorrelated wind
% velocity histories coherence model.

%
% Inputs:
%   Y   - Matrix of y-coordinates (lateral) of grid nodes. [Nzz,Nyy] = size(Y).
%   Z   - Matrix of z-coordinates (vertical) of grid nodes. Size must match Y.
%   U   - Vector or matrix of mean wind speed at each grid node. If a vector,
%         it must have one element for each node.
%   V   - Vector or matrix of mean wind speed at each grid node. If a vector,
%         it must have one element for each node.
%   W   - Vector or matrix of mean wind speed at each grid node. If a vector,
%         it must have one element for each node.
%   f   - Vector of frequency steps for the simulation.
%   u   - Matrix of uncorrelated wind velocities. Each column corresponds to
%         a grid node and each row to a time step.
%
% Outputs:
%   u_corr - Matrix of correlated wind velocities. Same dimensions as input
%            'u', with each column representing a grid node and each row
%            representing a time step.
%   v_corr - Matrix of correlated wind velocities. Same dimensions as input
%            'v', with each column representing a grid node and each row
%            representing a time step.
%   w_corr - Matrix of correlated wind velocities. Same dimensions as input
%            'w', with each column representing a grid node and each row
%            representing a time step.
%
% Original Author: E. Cheynet - uiB - Last modified 2024-04-16
% See also FFT, IFFT, CHOL, LDL.

% ##### Modified version ##### %

%% Prepare data for the coherence
% [Nzz,Nyy]=size(Y);
dy = abs(Y(:)'-Y(:)); % Matrix distance along y
dz = abs(Z(:)'-Z(:)); % Matrix distance along z
meanU = 0.5*(U(:)'+U(:)); % Mean wind velocity between each nodes

allthezs=Z(:).*ones(height(Z(:)),height(Y(:)));

%% Introduce coherence
% Define a frequency vector
f2s = [0 f(1:end-1)];
f2s = [f2s fliplr(f2s)];% 2-sided frequency array
tic
N = numel(f2s);

fftU = (fft(u));
fftV = (fft(v));
fftW = (fft(w));
cohU=NaN(size(allthezs));
cohV=NaN(size(allthezs));
% cohW=NaN(size(allthezs));

for ii=2:N
    for hh=1:height(allthezs)
        for bb=1:height(allthezs)
            [cohU(hh,bb)] = coh_OWN_U(meanU(hh,bb),dy(hh,bb),dz(hh,bb),f2s(ii),allthezs(hh,bb));
            [cohV(hh,bb)] = coh_OWN_V(meanU(hh,bb),dy(hh,bb),dz(hh,bb),f2s(ii),allthezs(hh,bb));
        end
    end
    if min(cohU(:))<0.98
        try
            C = chol(cohU,'lower');
            Cv = chol(cohV,'lower');
            % Cw = chol(cohW,'lower');
        catch exception
            [L,D]=ldl(cohU,'lower'); % a LDL decomposition is applied this time
            C = (L*sqrt(D));

            [Lv,Dv]=ldl(cohV,'lower'); % a LDL decomposition is applied this time
            Cv = (Lv*sqrt(Dv));

            % [Lw,Dw]=ldl(cohW,'lower'); % a LDL decomposition is applied this time
            % Cw = (Lw*sqrt(Dw));
        end
        fftU(ii,:)= C*exp(1i*angle(fftU(ii,:)')).*abs(fftU(ii,:)');

        fftV(ii,:)= Cv*exp(1i*angle(fftV(ii,:)')).*abs(fftV(ii,:)');

        % fftW(ii,:)= Cw*exp(1i*angle(fftW(ii,:)')).*abs(fftW(ii,:)');
    else
        fftU(ii,:)= mean(fftU(ii,:));
        fftV(ii,:)= mean(fftV(ii,:));
        fftW(ii,:)= mean(fftW(ii,:));
    end
end
toc
% zero mean value
fftU(1,:)=0;
fftV(1,:)=0;
fftW(1,:)=0;
% Get back to the time domain
u_corr = real(ifft(fftU));
v_corr = real(ifft(fftV));
w_corr = real(ifft(fftW));
% Add the mean wind speed
u_corr = (u_corr'+ U(:))';
v_corr = (v_corr'+ V(:))';
w_corr = (w_corr'+ W(:))';
    function [coh] = coh_OWN_U(meanU,dy,dz,f,z)

        % lateral separation
        ay_U     = 35;
        beta_y_U = 1.16;
        % vertical separation
        az_U     = 35;
        beta_z_U = 1.16;
        % Combine them into the coherence matrix for lateral and vertical
        % separations
        coh = exp(-sqrt((ay_U.*((f.*dy./meanU).^beta_y_U + (0.05.*dy./340.2).^beta_y_U)).^2 ...
          + (az_U.*((f.*dz./meanU).^beta_z_U + (0.05.*dz./340.2).^beta_z_U)).^2));
    end

    function [coh] = coh_OWN_V(meanU,dy,dz,f,z)

        % lateral separation
        ay_V      = 9;
        beta_y_V  = 1;
        % vertical separation
        az_V      = 9;
        beta_z_V  = 1;

        % Combine them into the coherence matrix for lateral and vertical
        % separations
        %coh = exp(-sqrt((ay_V.*(f.*dy./meanU).^beta_y_V).^2+(az_V.*(f.*dz./meanU).^beta_z_V).^2));
       coh = exp(-sqrt((ay_V.*((f.*dy./meanU).^beta_y_V + (0.05.*dy./340.2).^beta_y_V)).^2 ...
          + (az_V.*((f.*dz./meanU).^beta_z_V + (0.05.*dz./340.2).^beta_z_V)).^2));

    end
end
