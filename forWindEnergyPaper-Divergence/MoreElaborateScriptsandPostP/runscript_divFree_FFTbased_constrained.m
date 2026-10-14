clear
clc
close all
TOPFOLDERS = {'Data/KaimalDecay','Data/DiffDecay'};
TOPFOLDERS_DIV_FREE = {'Data/KaimalDecay_divfree_constrained','Data/DiffDecay_divfree_constrained'};
SUBFOLDERS = {'R_NABL'};
SEEDING = {'Seed1'};
tic;
for jj = 1:length(TOPFOLDERS)
    for ii = 1:length(SUBFOLDERS)
        for hh = 1:length(SEEDING)
            filepathhere = fullfile(TOPFOLDERS{jj}, SUBFOLDERS{ii}, SEEDING{hh}, 'forum_code_wnd_v2.wnd');
            filepathhere_DIV_FREE = fullfile(TOPFOLDERS_DIV_FREE{jj}, SUBFOLDERS{ii}, SEEDING{hh}, 'forum_code_wnd_v2.wnd');

            if ~exist(fileparts(filepathhere_DIV_FREE), 'dir') 
                mkdir(fileparts(filepathhere_DIV_FREE)); 
            end

            filepathhere_DIV_FREE_sum = fullfile(TOPFOLDERS_DIV_FREE{jj}, SUBFOLDERS{ii}, SEEDING{hh}, 'forum_code_wnd_v2');
            [velocity2D, y, z, ~, ~, ~, ~, dt, ~, ~, Summvars] = readfile_WND(filepathhere);
            filename = filepathhere_DIV_FREE;
            filenameSUM = filepathhere_DIV_FREE_sum;

            %%%%%%%%%% Correction begins %%%%%%%%%%%%%%%%%%%
            [Nt, ncomp, Ny, Nz] = size(velocity2D);
            U_mean = Summvars(3);   % Mean advection speed [m/s] 
            dx = U_mean * dt;       % Spatial step in x [Taylor's hypothesis]
            Nx = Nt;                % Let the number of x points equal the number of timesteps
            x = (0:Nx-1) * dx;      % x coordinate vector

            % Building a 3D velocity field (Nx, Ny, Nz, 3)
            velocity = buildVelocityField(velocity2D, Nx, Ny, Nz);

            % Grid spacings in y and z
            dy = mean(diff(y));
            dz = mean(diff(z));

            
            tolerance = 1e-9;
            max_iterations = 25;
            velocity_new = velocity;  % Initialize corrected velocity with the original
            div_v_og=computeDivergencefourth(velocity, Nx, Ny, Nz, dx, dy, dz);
            % keyboard
            for iter = 1:max_iterations
                % Compute divergence of current velocity field 
                div_v = computeDivergencefourth(velocity_new, Nx, Ny, Nz, dx, dy, dz);
                max_div = mean(abs(div_v(:)));
                fprintf('Iteration %d, mean divergence: %g\n', iter, max_div);

                if max_div < tolerance
                    break;
                end

                % --- FFT-based Poisson Solve and Spectral Differentiation ---
                % This function returns both the pressure (for diagnostics) and the gradient.
                [pressure, grad_phi] = solvePoissonAndGradientFFT(div_v, Nx, Ny, Nz, dx, dy, dz);

                % Update the velocity field to remove the divergence
                % Under relaxation factor of 0.9 used [one can play around with it if need be]
                velocity_new(:,:,:,1) = velocity_new(:,:,:,1) - 0.9 .* grad_phi(:,:,:,1);
                velocity_new(:,:,:,2) = velocity_new(:,:,:,2) - 0.9 .* grad_phi(:,:,:,2);
                velocity_new(:,:,:,3) = velocity_new(:,:,:,3) - 0.9 .* grad_phi(:,:,:,3);
                
                diff_u = velocity_new(:,:,:,1) - velocity(:,:,:,1);
                diff_v = velocity_new(:,:,:,2) - velocity(:,:,:,2);
                diff_w = velocity_new(:,:,:,3) - velocity(:,:,:,3);

                rmse_threshold_u = 0.35; 
                rmse_threshold_v = 0.35;
                rmse_threshold_w = 0.35;


                scaling_u = min(1, rmse_threshold_u ./ (abs(diff_u) + 1e-6)); % Avoid divide by zero
                velocity_new(:,:,:,1) = velocity(:,:,:,1) + scaling_u .* diff_u;

                scaling_v = min(1, rmse_threshold_v ./ (abs(diff_v) + 1e-6)); % Avoid divide by zero
                velocity_new(:,:,:,2) = velocity(:,:,:,2) + scaling_v .* diff_v;

                scaling_w = min(1, rmse_threshold_w ./ (abs(diff_w) + 1e-6)); % Avoid divide by zero
                velocity_new(:,:,:,3) = velocity(:,:,:,3) + scaling_w .* diff_w;
                

                
                iterdiv(iter) = mean(abs(div_v(:)));

                
            end
            toc;
            corrected_velocity(:,1,:,:) = real(velocity_new(:,:,:,1));
            corrected_velocity(:,2,:,:) = real(velocity_new(:,:,:,2));
            corrected_velocity(:,3,:,:) = real(velocity_new(:,:,:,3));

            % Please change below inputs as per your wind fields
            % Currently they are hard-coded to example wind fields supplied
            zHub = 170;
            z1 = 0;
            z = z - 5;
            
            numY = Ny;
            numZ = Nz;
            gridHt = 310;
            zOffset = z1 + gridHt/2;
            z0 = 0.0001;
            Clockwise = 1;
            UBAR = mean(corrected_velocity(:,1,18,18));
            TI_u = 100*(std(corrected_velocity(:,1,18,18)) / UBAR);
            TI_v = 100*(std(corrected_velocity(:,2,18,18)) / UBAR);
            TI_w = 100*(std(corrected_velocity(:,3,18,18)) / UBAR);
            SummVars = [zHub; Clockwise; UBAR; TI_u; TI_v; TI_w];
            
            writeBLgrid(filename, corrected_velocity, dy, dz, dt, zOffset, z0, SummVars)
            write_wind_file(dt, UBAR, TI_u, TI_v, TI_w, filenameSUM)
        end
    end
end
%% Quick post-processing example (Sanity check)
% histogram(div_v_og(:),Normalization="probability",BinWidth=0.01)
% hold on
% histogram(div_v(:),Normalization="probability",BinWidth=0.01)



%% ==========================
%% Function Definitions
%% ==========================

% Function to build a 3D velocity field from 2D data
function velocity = buildVelocityField(velocity2D, Nx, Ny, Nz)
    velocity = zeros(Nx, Ny, Nz, 3);
    velocity(:,:,:,1) = velocity2D(:,1,:,:);  % u component (x-direction)
    velocity(:,:,:,2) = velocity2D(:,2,:,:);  % v component (y-direction)
    velocity(:,:,:,3) = velocity2D(:,3,:,:);  % w component (z-direction)
end

% Compute divergence using fourth-order centered differences
function div_v = computeDivergencefourth(velocity, Nx, Ny, Nz, dx, dy, dz)
    div_v = zeros(Nx, Ny, Nz);
    for k = 1:Nz
        for j = 1:Ny
            for i = 1:Nx
                % du/dx (Fourth-order centered difference)
                if i <= 2 || i >= Nx-1
                    du_dx = 0; % Boundary condition for first and last points
                else
                    du_dx = (-velocity(i+2, j, k, 1) + 8*velocity(i+1, j, k, 1) - 8*velocity(i-1, j, k, 1) + velocity(i-2, j, k, 1)) / (12 * dx);
                end

                % dv/dy (Fourth-order periodic difference)
                jp = mod(j, Ny) + 1;
                jm = mod(j-2, Ny) + 1;
                jp2 = mod(j+1, Ny) + 1;
                jm2 = mod(j-3, Ny) + 1;
                dv_dy = (-velocity(i, jp2, k, 2) + 8*velocity(i, jp, k, 2) - 8*velocity(i, jm, k, 2) + velocity(i, jm2, k, 2)) / (12 * dy);

                % dw/dz (Fourth-order with no-flow BC)
                if k <= 2 || k >= Nz-1
                    dw_dz = 0;
                else
                    dw_dz = (-velocity(i, j, k+2, 3) + 8*velocity(i, j, k+1, 3) - 8*velocity(i, j, k-1, 3) + velocity(i, j, k-2, 3)) / (12 * dz);
                end

                div_v(i, j, k) = du_dx + dv_dy + dw_dz;
            end
        end
    end
end

%-------------------------------------------------------
% FFT-based Poisson solver and spectral gradient computation
%-------------------------------------------------------
function [pressure, grad_phi] = solvePoissonAndGradientFFT(div_v, Nx, Ny, Nz, dx, dy, dz)
    % Step 1: Extend div_v for Neumann BCs by mirroring in x and z
    div_v_extended = [div_v; flip(div_v, 1)];  % Mirror in x-direction
    div_v_extended = cat(3, div_v_extended, flip(div_v_extended, 3));  % Mirror in z-direction

    % Extended domain sizes
    Nx_ext = 2 * Nx;
    Nz_ext = 2 * Nz;

    % Step 2: Compute FFT of the extended divergence field
    div_v_hat = fftn(div_v_extended, [Nx_ext, Ny, Nz_ext]);

    % Step 3: Define wavenumber vectors for the extended grid
    kx_vec = (2 * pi / (Nx_ext * dx)) * ifftshift(-floor(Nx_ext/2) : ceil(Nx_ext/2)-1);
    ky_vec = (2 * pi / (Ny * dy)) * ifftshift(-floor(Ny/2) : ceil(Ny/2)-1);
    kz_vec = (2 * pi / (Nz_ext * dz)) * ifftshift(-floor(Nz_ext/2) : ceil(Nz_ext/2)-1);

    [kx, ky, kz] = ndgrid(kx_vec, ky_vec, kz_vec);

    % Step 4: Compute squared wavenumbers (avoid division by zero)
    k2 = kx.^2 + ky.^2 + kz.^2;
    k2(1,1,1) = 1;  % To prevent division by zero for the zero-frequency (DC) mode

    % Step 5: Solve Poisson equation in Fourier space
    pressure_hat = -div_v_hat ./ k2;
    pressure_hat(1,1,1) = 0;  % Enforce zero-mean (remove arbitrary constant)

    % Step 6: Inverse FFT to obtain pressure in the extended domain
    pressure_extended = real(ifftn(pressure_hat));

    % Step 7: Extract the original pressure field
    pressure = pressure_extended(1:Nx, :, 1:Nz);

    % Step 8: Compute the gradient in Fourier space
    grad_phi_x_hat = 1i * kx .* pressure_hat;
    grad_phi_y_hat = 1i * ky .* pressure_hat;
    grad_phi_z_hat = 1i * kz .* pressure_hat;

    % Step 9: Inverse FFT to get the gradient in real space (extended domain)
    grad_phi_x_extended = real(ifftn(grad_phi_x_hat));
    grad_phi_y_extended = real(ifftn(grad_phi_y_hat));
    grad_phi_z_extended = real(ifftn(grad_phi_z_hat));

    % Step 10: Extract the gradients for the original domain
    grad_phi_x = grad_phi_x_extended(1:Nx, :, 1:Nz);
    grad_phi_y = grad_phi_y_extended(1:Nx, :, 1:Nz);
    grad_phi_z = grad_phi_z_extended(1:Nx, :, 1:Nz);

    % Combine the three components into a 4D array (Nx, Ny, Nz, 3)
    grad_phi = cat(4, grad_phi_x, grad_phi_y, grad_phi_z);

end

%%
function [velocity, y, z, nz, ny, dz, dy, dt, zHub, z1, SummVars] = readfile_WND(FileName)
%[velocity, y, z, nz, ny, dz, dy, dt, zHub, z1, SummVars] = readfile_WND(FileName)
% formerly readBLgrid()
% Input:
% FileName       - string, containing file name to open (.wnd extension is optional)
%
% Output:
%  velocity      - 4-D vector: time, velocity component, iy, iz
%  y             - 1-D vector: horizontal locations y(iy)
%  z             - 1-D vector: vertical locations z(iz)
%  nz, ny        - scalars: number of points in the vertical/horizontal
%                  direction of the grid
%  dz, dy, dt    - scalars: distance between two points in the vertical [m]/
%                  horizontal [m]/time [s] dimension
% zHub           - hub height [m]
% z1             - vertical location of bottom of grid [m above ground level]
% SumVars        - variables from the summary file (zHub, Clockwise, UBAR, TI_u, TI_v, TI_w)

%%
len    = length(FileName);
ending = FileName(len-3:len);

if strcmpi( ending, '.wnd' )
    FileName = FileName(1:len-4);
end

%-------------------------------------------------------------

% initialize variables
fileFmt  = 'int16';
ConvFact = 1.0; %results in meters and seconds

str      = {'HUB HEIGHT','CLOCKWISE','UBAR','TI(U','TI(V','TI(W'};  %MUST be in UPPER case
numVars  = length(str);
SummVars = zeros(numVars, 1);

%% -----------------------------------------
%  READ THE HEADER OF THE BINARY FILE
%  -----------------------------------------
fid_wnd   = fopen( [ FileName '.wnd' ] );
if ( fid_wnd <= 0 )
    error( 'Wind file could not be opened.' );
end

nffc  = fread( fid_wnd, 1, 'int16' );                     % number of components

if nffc ~= -99  % AN OLD-STYLE AERODYN WIND FILE
    dz      = fread( fid_wnd, 1, 'int16' );               % delta z in mm
    dy      = fread( fid_wnd, 1, 'int16' );               % delta y in mm
    dx      = fread( fid_wnd, 1, 'int16' );               % delta x (actually t in this case) in mm
    nt      = fread( fid_wnd, 1, 'int16' );               % half number of time steps
    MFFWS   = fread( fid_wnd, 1, 'int16' );               % 10 times mean FF wind speed, should be equal to MWS
    fread( fid_wnd, 5, 'int16' );               % unnecessary lines
    nz      = fread( fid_wnd, 1, 'int16' );               % 1000 times number of points in vertical direction, max 32
    ny      = fread( fid_wnd, 1, 'int16' );               % 1000 times the number of points in horizontal direction, max 32
    fread( fid_wnd, 3*(-nffc-1), 'int16' );

    % convert the integers to real numbers
    nffc     = -nffc;
    dz       = 0.001*ConvFact*dz;
    dy       = 0.001*ConvFact*dy;
    dx       = 0.001*ConvFact*dx;
    MFFWS    = 0.1*ConvFact*MFFWS;
    nz       = fix( mod(nz,2^16) / 1000 );                % the mod 2^16 is a work around for somewhat larger grids
    ny       = fix( mod(ny,2^16) / 1000 );                % the mod 2^16 is a work around for somewhat larger grids

else %== -99, THE NEWER-STYLE AERODYN WIND FILE
    fc       = fread( fid_wnd, 1, 'int16' );              % should be 4 to allow turbulence intensity to be stored in the header

    if fc == 4
        nffc     = fread( fid_wnd, 1, 'int32' );              % number of components (should be 3)
        lat      = fread( fid_wnd, 1, 'float32' );            % latitude (deg)
        z0       = fread( fid_wnd, 1, 'float32' );            % Roughness length (m)
        zOffset  = fread( fid_wnd, 1, 'float32' );            % Reference height (m) = Z(1) + GridHeight / 2.0
        TI_U     = fread( fid_wnd, 1, 'float32' );            % Turbulence Intensity of u component (%)
        TI_V     = fread( fid_wnd, 1, 'float32' );            % Turbulence Intensity of v component (%)
        TI_W     = fread( fid_wnd, 1, 'float32' );            % Turbulence Intensity of w component (%)
    else
        if fc > 2
            nffc = 3;
        else
            nffc = 1;
        end
        TI_U  = 1;
        TI_V  = 1;
        TI_W  = 1;

        if fc == 8 ... %MANN model
                || fc == 7 % General Kaimal
            HeadRec = fread(fid_wnd,1,'int32');
            nffc    = fread(fid_wnd,1,'int32');  %nffc?
        end

    end %fc == 4
    %%
    dz       = fread( fid_wnd, 1, 'float32' );            % delta z in m
    dy       = fread( fid_wnd, 1, 'float32' );            % delta y in m
    dx       = fread( fid_wnd, 1, 'float32' );            % delta x in m
    nt       = fread( fid_wnd, 1, 'int32' );              % half the number of time steps
    MFFWS    = fread( fid_wnd, 1, 'float32');             % mean full-field wind speed

    fread( fid_wnd, 3, 'float32' );            % zLu, yLu, xLu: unused variables (for BLADED)
    fread( fid_wnd, 2, 'int32' );              % unused variables (for BLADED) [unused integer, random seed]
    nz       = fread( fid_wnd, 1, 'int32' );              % number of points in vertical direction
    ny       = fread( fid_wnd, 1, 'int32' );              % number of points in horizontal direction
    if (nffc==3)
        fread( fid_wnd, 2*nffc, 'int32' );     % other length scales: unused variables (for BLADED)
    end

    if fc == 7
        CohDec = fread(fid_wnd,1,'float32');
        CohLc  = fread(fid_wnd,1,'float32');
    elseif fc == 8        % MANN model
        gamma  = fread(fid_wnd,1,'float32');               % MANN model shear parameter
        Scale  = fread(fid_wnd,1,'float32');               % MANN model scale length
        fread(fid_wnd,4,'float32');
        fread(fid_wnd,3,'int32');
        fread(fid_wnd,2,'float32');
        fread(fid_wnd,3,'int32');
        fread(fid_wnd,2,'float32');
    end

    SummVars(3:6) = [MFFWS, TI_U, TI_V, TI_W];

end % old or new bladed styles

nt     = max([nt*2,1]);
dt     = dx/MFFWS;

%% -----------------------------------------
%  READ THE SUMMARY FILE FOR SCALING FACTORS
%  -----------------------------------------
disp('Reading the summary file....');

indx     = SummVars;
fid_sum  = fopen( [ FileName '.sum' ] );

if ( fid_sum <= 0 )
    fclose(fid_sum);
    fclose(fid_wnd);

    error(['Could not open the summary file: ' FileName '.sum']);
end

while ( any( indx == 0 ) )  %MFFWS and the TIs should not be zero
    line  = fgetl(fid_sum);
    if ~ischar(line)
        % We reached the end of the file

        fclose(fid_sum);
        fclose(fid_wnd);

        error('Reached the end of summary file without all necessary data.');
    end

    line  = upper(line);
    findx = strfind(line,'=')+1;        %first index
    if isempty(findx)
        findx = 1;
    end
    lindx = length(line);               %last  index

    i = 1;
    while i <= numVars
        if indx(i)==0

            k = findstr(line, str{i} );
            if ~isempty(k)              % we found a string we're looking for
                indx(i) = k;
                k=strfind(line,'%');
                if ~isempty(k)
                    lindx = max(findx,k-1);
                end

                tmp = strtok(line(findx:lindx));
                if isempty( str2num(tmp) ) %bjj: was str2double, but that returns NaN for 'T'/'F'
                    if strcmpi( tmp(1), 'T' )
                        SummVars(i) = 1;
                    else
                        SummVars(i) = -1;  %use this for false instead of zero.
                    end
                else
                    SummVars(i) = str2double( tmp );
                    break;
                end %if isempty( str2num(tmp) )
            end % ~isempty(k)
        end % indx(i)==0
        i = i + 1;
    end %while i

end %while any(indx==0)

% read the rest of the file to get the grid height offset, if it's there
ZGoffset = 0.0;

while ( true )
    line  = fgetl(fid_sum);

    if ~ischar(line)
        break;
    end

    line  = upper(line);
    findx = strfind(line,'HEIGHT OFFSET');

    if ~isempty(findx)
        lindx = length(line);
        findx = strfind(line,'=')+1;
        if isempty(findx)
            findx = 1;
        end
        ZGoffset = str2double( strtok(line(findx:lindx))); %z grid offset
        break;
    end
end

LHR = false; % left-hand rule (to flip sign for v wind component)
while ( true )
    line  = fgetl(fid_sum);

    if ~ischar(line)
        break;
    end

    line  = upper(line);
    findx = strfind(line,'BLADED LEFT-HAND RULE');

    if ~isempty(findx)
        LHR = true;
        break;
    end
end

fclose(fid_sum);
%%
%-----------------------------------------
%READ THE GRID DATA FROM THE BINARY FILE
%-----------------------------------------
disp('Reading and scaling the grid data...');

% nffc     = 3;
nv       = nffc*ny*nz;               % the size of one time step
Scale    = 0.00001*SummVars(3)*SummVars(4:6);
Offset   = [SummVars(3) 0 0];
if LHR % Bladed defined the v-component opposite of the right-hand rule...
    Scale(2) = -Scale(2);
end

%%
velocity = zeros(nt,nffc,ny,nz);

%bjj: disp( [nv nffc ny nz nt] );

if SummVars(2) > 0 %clockwise rotation
    %flip the y direction....
    %let's change the dimension of velocity so that it's 4-d instead of 3-d
    y_ix = ny:-1:1;
else
    y_ix = 1:ny;
end

% [v cnt] = fread( fid_wnd, nv, fileFmt );
% if cnt < nv
%     error(['Could not read entire file: at grid record ' num2str( (it-1)*nv+cnt2 ) ' of ' num2str(nrecs)]);
% end
% disp('Scaling the grid data...');

for it = 1:nt

    [v cnt] = fread( fid_wnd, nv, fileFmt );
    if cnt < nv
        error(['Could not read entire file: at grid record ' num2str( (it-1)*nv+cnt ) ' of ' num2str(nv*nt)]);
    end

    cnt2 = 1;
    for iz = 1:nz
        for iy = y_ix
            for k=1:nffc
                velocity(it,k,iy,iz) = v(cnt2)*Scale(k) + Offset(k);
                cnt2 = cnt2 + 1;
            end %for k
        end %iy
    end % iz

end %it

fclose(fid_wnd);


y    = (0:ny-1)*dy - dy*(ny-1)/2;

zHub = SummVars(1);
z1   = zHub - ZGoffset - dz*(nz-1)/2;  %this is the bottom of the grid
z    = (0:nz-1)*dz + z1;

disp('Finished.');
disp('');

return;
end
%%
function writeBLgrid(FileName, velocity, dy, dz, dt, zOffset, z0, SummVars)


% writeBLgrid('testWind', velocity, y, z, nz, ny, dz, dy, dt, zHub, z1, SummVars)
%
% writes wind velocity data to binary .wnd files
%

% inputs - FileName: Name of output .wnd-file (extension will be added)
%        - velocity: 4D-array: (time, 3D-windcomp, y, z)
%        - dz
%        - dy
%        - dt
%        - SumVars: 6 variables from the summary file {zHub, Clockwise, UBAR, TI_u, TI_v, TI_w}
%                                                        90          1    12    12   9.6     6
%        - zOffset: (=90) Reference height (m) = Z(1) + GridHeight / 2.0
%        - z0 = .03;      Roughness length (m)

%-----------------------------------------
% INITIALIZE VARIABLES
%-----------------------------------------

fileFmt  = 'int16';

fc = 4;                 % should be 4 to allow turbulence intensity to be stored in the header
lat = 0;                % latitude (deg)

zHub = SummVars(1);     % hub height [m]
MFFWS = SummVars(3);    % mean full-field wind speed
TI_U = SummVars(4);     % Turbulence Intensity of u component (%)
TI_V = SummVars(5);     % Turbulence Intensity of v component (%)
TI_W = SummVars(6);     % Turbulence Intensity of w component (%)

[nt, nffc, ny, nz] = size(velocity); % determin dimensions of windfield, e.g. [1286,3,23,23]

z1 =5;  %this is the bottom of the grid

dx = dt*MFFWS;          % delta x in m
nt_header = floor(nt/2);% half the number of time steps


%-----------------------------------------
% OPEN FILE
%-----------------------------------------
len    = length(FileName);
ending = FileName(len-3:len);

if strcmpi( ending, '.wnd' )
    FileName = FileName(1:len-4);
end
fid_wnd   = fopen( [ FileName '.wnd' ], 'w' ); % Open file, or create new file, for writing; discard existing contents, if any.
if ( fid_wnd <= 0 )
    error( 'Wind file could not be opened.' );
    return;
end


%-----------------------------------------
%WRITE THE HEADER OF THE BINARY FILE
%-----------------------------------------
% THE NEWER-STYLE AERODYN WIND FILE

fwrite( fid_wnd, -99, 'int16' );             % number of components
fwrite( fid_wnd, fc, 'int16' );              % should be 4 to allow turbulence intensity to be stored in the header
fwrite( fid_wnd, nffc, 'int32' );            % number of components (should be 3)
fwrite( fid_wnd, lat, 'float32' );           % latitude (deg)
fwrite( fid_wnd, z0, 'float32' );            % Roughness length (m)
fwrite( fid_wnd, zOffset, 'float32' );       % Reference height (m) = Z(1) + GridHeight / 2.0
fwrite( fid_wnd, TI_U, 'float32' );          % Turbulence Intensity of u component (%)
fwrite( fid_wnd, TI_V, 'float32' );          % Turbulence Intensity of v component (%)
fwrite( fid_wnd, TI_W, 'float32' );          % Turbulence Intensity of w component (%)

fwrite( fid_wnd, dz, 'float32' );            % delta z in m
fwrite( fid_wnd, dy, 'float32' );            % delta y in m
fwrite( fid_wnd, dx, 'float32' );            % delta x in m
fwrite( fid_wnd, nt_header, 'int32' );       % half the number of time steps
fwrite( fid_wnd, MFFWS, 'float32');          % mean full-field wind speed

fwrite( fid_wnd, [0 0 0], 'float32' );    % unused variables (for BLADED): write zeros
fwrite( fid_wnd, [0 0], 'int32' );        % unused variables (for BLADED): write zeros
fwrite( fid_wnd, nz, 'int32' );              % number of points in vertical direction
fwrite( fid_wnd, ny, 'int32' );              % number of points in horizontal direction
fwrite( fid_wnd, zeros(3*(nffc-1),1), 'int32' );     % unused variables (for BLADED): write zeros
%-----------------------------------------
% WRITE GRID DATA
%-----------------------------------------
Scale    = 0.00001*SummVars(3)*SummVars(4:6);
Offset   = [SummVars(3) 0 0];


if SummVars(2) > 0 %clockwise rotation
    %flip the y direction....
    %let's change the dimension of velocity so that it's 4-d instead of 3-d
    y_ix = ny:-1:1;
else
    y_ix = 1:ny;
end

v = zeros(nz*ny*nffc,1);

for it = 1:nt
    cnt = 1;
    for iz = 1:nz
        for iy = y_ix
            for k=1:nffc
                v(cnt) = (velocity(it,k,iy,iz) - Offset(k))/Scale(k);
                cnt = cnt + 1;
                %velocity(it,k,iy,iz) = v(cnt2)*Scale(k) + Offset(k);
            end %for k
        end %iy
    end % iz

    fwrite( fid_wnd, v, fileFmt );

end %it

%-----------------------------------------
% CLOSE .WND FILE
%-----------------------------------------
fclose(fid_wnd);



% y    = [0:ny-1]*dy - dy*(ny-1)/2;
% zHub = SummVars(1);
% z1   = zHub - ZGoffset - dz*(nz-1)/2;  %this is the bottom of the grid
% z    = [0:nz-1]*dz + z1;
end
%%
function write_wind_file(dt, UBar, TI_u, TI_v, TI_w,filename)
% Open file for writing
% fileID = fopen('forum_code_wnd_v2.sum', 'w');

% pathsaving = ['./', char(subfolder)];


fileID = fopen([filename,'.sum'], 'w');


% Write header
fprintf(fileID, 'T  Clockwise\n');
fprintf(fileID, '0 -  NONE  IEC turbulence models scaled to exact specified standard deviation\n\n');
fprintf(fileID, 'Turbine/Model Specifications:\n');
fprintf(fileID, '         35  Vertical grid-point matrix dimension\n');
fprintf(fileID, '         35  Horizontal grid-point matrix dimension\n');
fprintf(fileID, '         %.2f  Time step [seconds]\n', dt);
fprintf(fileID, '         700  Analysis time [seconds]\n');
fprintf(fileID, '         700  Usable output time [seconds]\n\n');

% Write hub height
fprintf(fileID, '170.000  Hub height [m]\n\n');

% Write UBar and TI values
fprintf(fileID, 'UBar  =  %.4f m/s\n', UBar);
fprintf(fileID, 'TI(u) =  %.4f %%\n', TI_u);
fprintf(fileID, 'TI(v) =  %.4f %%\n', TI_v);
fprintf(fileID, 'TI(w) =  %.4f %%\n', TI_w);

% Close file
fclose(fileID);
end
