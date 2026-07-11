function [u_corr,v_corr,w_corr] = func_windSimFaster_IECKAI(Y,Z,U,V,W,f,u,v,w,vhubplease)
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
%% Prepare data for the coherence
% [Nzz,Nyy]=size(Y);

vhub=vhubplease;

dy = abs(Y(:)'-Y(:)); % Matrix distance along y
dz = abs(Z(:)'-Z(:)); % Matrix distance along z

% Combine them into the coherence matrix for lateral and vertical
% separations
K = -sqrt(dy.^2+dz.^2);
% Anonymous function
cohIECKai = @(K,f) (exp(-12.*sqrt((f.*K./vhub).^2+(0.12.*K./340.2).^2))).^2;


%% Introduce coherence
% Define a frequency vector
f2s = [0 f(1:end-1)];
f2s = [f2s fliplr(f2s)];% 2-sided frequency array
tic
N = numel(f2s);

fftU = (fft(u));
fftV = (fft(v));
fftW = (fft(w));
for ii=2:N
    cohU = cohIECKai(K,f2s(ii));
    if min(cohU(:))<0.98
        try
            C = chol(cohU,'lower');
        catch exception
            [L,D]=ldl(cohU,'lower'); % a LDL decomposition is applied this time
            C = (L*sqrt(D));
        end
        fftU(ii,:)= C*exp(1i*angle(fftU(ii,:)')).*abs(fftU(ii,:)');
        fftV(ii,:)= C*exp(1i*angle(fftV(ii,:)')).*abs(fftV(ii,:)');
        % fftW(ii,:)= C*exp(1i*angle(fftW(ii,:)')).*abs(fftW(ii,:)');
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
end