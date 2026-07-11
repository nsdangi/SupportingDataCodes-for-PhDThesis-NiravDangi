clear
clc
addpath("functions");
% addpath("code_syntheticwf\CUSTOM\functions")

rng_seed=123; % Seed number set

zhere        = linspace(40,240,21); % example altitudes
z0           = 0.01;  % example surface roughness
hubheight    = 140;   % example hub height or z_ref
u_ref        = 10;    % example hub height u component
v_ref        = 1 ;    % example hub height v component
w_ref        = 1 ;    % example hub height w component

TI_u_ref     = 7 ;    % example hub height TI u component
TI_v_ref     = 2 ;    % example hub height TI v component
TI_w_ref     = 2 ;    % example hub height TI w component
 

[~, hubheightIDX] = min(abs(zhere - hubheight));

% stabIDX = 1; % Very unstable (VU)
% stabIDX = 2; % Slightly unstable (SU)
% stabIDX = 3; % Neutral (N)
% stabIDX = 4; % Slightly stable (SS)
% stabIDX = 5; % Very stable (VS)

stabIDX = 5; 

powerlawexp = 0.25;

% Power law wind shear input 
% OR set known wind shear manually
uhere= func_powerlawshear(u_ref, hubheight, zhere, powerlawexp);

% Example LMO values
LMO_map = [-100, -350, 1e10, 350, 100]; % VU, SU, N, SS, VS
LMO = LMO_map(stabIDX);           % assign LMO based on stability index

% Example Inversion height values
ZIin_map = [1500, 700, 1000, 300, 200]; % VU, SU, N, SS, VS
ZIin  = ZIin_map(stabIDX);  

aslheight = 150; % Fitting M-O profile upto what altitudes
[~, aslheightIDX] = min(abs(zhere - aslheight));


% Find u star from shear and M-O length
params = func_fit_MO_profile_knownL(zhere(1:aslheightIDX), uhere(1:aslheightIDX), z0, LMO);
ustar=params.ustar;

% Wind shear with M-O formulations
fittedu = func_wind_profile_MO(zhere, z0, ustar, LMO);

ucomp=fittedu;
vcomp=ones(size(fittedu)).*v_ref;
wcomp=ones(size(fittedu)).*w_ref;

vhubplease=fittedu(hubheightIDX);

Nyy=length(zhere);
Nzz=length(zhere);

gridHt=max(zhere)-min(zhere);

dt = 0.5;

% Current implementation generates time-series of length 2^(timeexp)
% 2^10 --> 1024 time-steps -->if dt= 0.5 --> total time 512 seconds

timeexp=10;



opfolderName = 'Outputs';
mkdir(opfolderName);


output_fname_bts = ['Outputs/example_owncoh.bts'];
output_fname_wnd = ['Outputs/example_owncoh'];

runWindSim_OWNcoh(hubheight, vhubplease, ucomp',vcomp', wcomp',zhere, ...
    TI_u_ref, TI_v_ref, TI_w_ref, Nyy, Nzz, gridHt, timeexp, z0, rng_seed, ...
    output_fname_bts, output_fname_wnd,stabIDX,ustar,LMO,ZIin,dt);



%%
function runWindSim_OWNcoh(hubheight, vhubplease, ucomp,vcomp, wcomp,heightsarr, ...
    ti_u, ti_v, ti_w, Nyy, Nzz, gridHt, timeexp, z0, rng_seed, ...
    output_fname_bts, output_fname_wnd,stabIDX,u_star,Lmo,ZIin,dt)
% runWindSim  - Wind turbulence simulation wrapper 
%
% Syntax:
%   runWindSim(hubheight, vhubplease, vcomp, wcomp, ...
%              ti_u, ti_v, ti_w, shear_exp, Nyy, Nzz, gridHt, timeexp, z0, rng_seed, ...
%              output_fname_bts, output_fname_wnd)
%
% Inputs:
%   hubheight       - Hub height [m]
%   vhubplease      - Mean wind speed at hub [m/s]
%   vcomp, wcomp    - Mean V, W components [m/s]
%   ti_u, ti_v, ti_w- Turbulence intensities [%]
%   shear_exp       - Shear exponent
%   Nyy, Nzz        - Grid size in y and z directions
%   gridHt          - Vertical extent of the grid [m]
%   timeexp         - Exponent for # of timesteps (N=2^timeexp)
%   z0              - Surface roughness length [m]
%   rng_seed        - Random number generator seed
%   output_fname_bts- Output filename for .bts file
%   output_fname_wnd- Output filename prefix for wind grid

%% --- Derived parameters ---
sig_u = ti_u*vhubplease/100;
sig_v = ti_v*vhubplease/100;
sig_w = ti_w*vhubplease/100;

kappa = 0.4; % von Karman constant
spany = linspace(-round(gridHt/2),round(gridHt/2),Nyy);
spanz = heightsarr;


rng(rng_seed)  % Set random seed

%% --- Grid and shear profile ---
[Y,Z] = meshgrid(spany, spanz);
U     = ones(size(Y)).*ucomp;
V     = ones(size(Y)).*vcomp;
W     = ones(size(Y)).*wcomp;

%% --- Time/frequency setup ---
fs = 1/dt;
[t,f] = func_getSamplingPara(timeexp,fs);
N = numel(t);
Nfreq = numel(f);

%% --- Spectra + random processes ---
Mgrid = numel(Z(:));
u = zeros(N,Mgrid);
v = zeros(N,Mgrid);
w = zeros(N,Mgrid);

Su = zeros(Nfreq,Mgrid);
Sv = zeros(Nfreq,Mgrid);
Sw = zeros(Nfreq,Mgrid);

dummyU = U(:);
dummyZ = Z(:);

for ii = 1:Mgrid
    

    if stabIDX ==1 || stabIDX==2
        [Su(:,ii),Sv(:,ii),Sw(:,ii),~,~,~] = ...
        func_KaimalModel_USABL(dummyU(ii),dummyZ(ii),f,u_star,vhubplease,sig_u,Lmo,ZIin);


    elseif stabIDX ==3
        [Su(:,ii),Sv(:,ii),Sw(:,ii),~,~,~] = ...
        func_KaimalModel_NABL(dummyU(ii),dummyZ(ii),f,u_star,vhubplease,sig_u);

    elseif stabIDX==4 || stabIDX==5

        [Su(:,ii),Sv(:,ii),Sw(:,ii),~,~,~] = ...
        func_KaimalModel_SABL(dummyU(ii),dummyZ(ii),f,u_star,vhubplease,sig_u,Lmo);

    end

    [u(:,ii),t] = func_randomProcess(f,Su(:,ii));
    [v(:,ii),t] = func_randomProcess(f,Sv(:,ii));
    [w(:,ii),t] = func_randomProcess(f,Sw(:,ii));
    
end
%% --- Correlated wind fields ---

if stabIDX ==1 || stabIDX==2
    [u_corr,v_corr,w_corr] = func_windSimFaster_USABL(Y,Z,U,V,W,f,u,v,w);


elseif stabIDX ==3
    [u_corr,v_corr,w_corr] = func_windSimFaster_NABL(Y,Z,U,V,W,f,u,v,w);

elseif stabIDX==4 || stabIDX==5

    [u_corr,v_corr,w_corr] = func_windSimFaster_SABL(Y,Z,U,V,W,f,u,v,w);
    


end
u_corr = reshape(u_corr,[N,Nzz,Nyy]);
v_corr = reshape(v_corr,[N,Nzz,Nyy]);
w_corr = reshape(w_corr,[N,Nzz,Nyy]);

%% --- Hub location indices based on physical coordinates ---
y = Y(1,:)';      % lateral grid
z = Z(:,1);       % vertical grid
[~, idx_y] = min(abs(y - 0));          % hub at y=0
[~, idx_z] = min(abs(z - hubheight));  % hub height

%% --- Extract hub velocity using indices ---
uHub = u_corr(:,idx_z,idx_y);
timemax = size(u_corr,1)*dt;

%% --- Write BTS output ---
ntower = 0;
func_write2bts(output_fname_bts,u_corr,v_corr,w_corr,t,y,z,hubheight,mean(uHub),ntower);

%% --- Post-process: read back ---
[velocity1, ~, y, z, ~, nz, ny, dz, dy, dt, hubheight, z1,mffws] = ...
    func_readfile_BTS(output_fname_bts);

sTD_hub_diffU = std(squeeze(velocity1(:, 1, idx_y,idx_z)));  % For U component
sTD_hub_diffV = std(squeeze(velocity1(:, 2, idx_y,idx_z)));  % For V component
sTD_hub_diffW = std(squeeze(velocity1(:, 3, idx_y,idx_z)));  % For W component

% Calculate mean velocity field across time for U, V, and W components
meanVELDIFF_U = mean(velocity1(:, 1, :, :), 1);
meanVELDIFF_V = mean(velocity1(:, 2, :, :), 1);
meanVELDIFF_W = mean(velocity1(:, 3, :, :), 1);

sTD_hub_targetU=sig_u;
sTD_hub_targetV=sig_v;
sTD_hub_targetW=sig_w;

% Calculate scaling factors for U, V, and W components
scalingDIFF_U = sTD_hub_targetU / sTD_hub_diffU;
scalingDIFF_V = sTD_hub_targetV / sTD_hub_diffV;
scalingDIFF_W = sTD_hub_targetW / sTD_hub_diffW;

% Apply scaling to velocity fluctuations for each component
scaledVelocity_U = meanVELDIFF_U + (velocity1(:, 1, :, :) - meanVELDIFF_U) * scalingDIFF_U;
scaledVelocity_V = meanVELDIFF_V + (velocity1(:, 2, :, :) - meanVELDIFF_V) * scalingDIFF_V;
scaledVelocity_W = meanVELDIFF_W + (velocity1(:, 3, :, :) - meanVELDIFF_W) * scalingDIFF_W;

scaledVelocity = zeros(size(velocity1, 1), 3, size(velocity1, 3), size(velocity1, 4));
scaledVelocity(:, 1, :, :) = scaledVelocity_U;  % U component
scaledVelocity(:, 2, :, :) = scaledVelocity_V;  % V component
scaledVelocity(:, 3, :, :) = scaledVelocity_W;  % W component

meanvel = scaledVelocity;
zOffset = z1+gridHt/2;

%% --- TI & summary at hub location ---
i = idx_z;
j = idx_y;

Clockwise = 1;
UBAR = mffws;
TI_u = 100*(std(meanvel(:,1,i,j))/mean(meanvel(:,1,i,j)));
TI_v = 100*(std(meanvel(:,2,i,j))/mean(meanvel(:,1,i,j)));
TI_w = 100*(std(meanvel(:,3,i,j))/mean(meanvel(:,1,i,j)));
SummVars = [hubheight; Clockwise; UBAR; TI_u; TI_v; TI_w];

%% --- Write outputs ---
func_writeBLgrid(output_fname_wnd, meanvel, dy, dz, dt, zOffset, z0, SummVars)
func_write_wind_file(dt, UBAR, TI_u, TI_v, TI_w, output_fname_wnd,Nyy,Nzz,timemax,hubheight)

end
%%







